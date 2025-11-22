#include <algorithm>
#include <cctype>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using Byte = std::uint8_t;

enum class BerClass : std::uint8_t {
    Universal = 0,
    Application = 1,
    ContextSpecific = 2,
    Private = 3
};

struct BerTag {
    BerClass cls = BerClass::Universal;
    bool constructed = false;
    std::uint32_t number = 0;
};

class PinParseError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

std::vector<Byte> ReadFile(const std::filesystem::path &path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Unable to open file: " + path.string());
    }

    return std::vector<Byte>(std::istreambuf_iterator<char>(in),
                             std::istreambuf_iterator<char>());
}

std::uint32_t ReadBE32(const std::vector<Byte> &buffer, std::size_t &offset)
{
    if (offset + 4 > buffer.size()) {
        throw PinParseError("Unexpected end of file while reading 32-bit value");
    }

    std::uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        value = (value << 8) | buffer[offset + i];
    }
    offset += 4;
    return value;
}

BerTag ReadTag(const std::vector<Byte> &buffer, std::size_t &offset)
{
    if (offset >= buffer.size()) {
        throw PinParseError("Unexpected end of buffer while reading BER tag");
    }

    const Byte first = buffer[offset++];
    BerTag tag;
    tag.cls = static_cast<BerClass>((first & 0b1100'0000) >> 6);
    tag.constructed = (first & 0b0010'0000) != 0;
    tag.number = first & 0b0001'1111;

    if (tag.number == 0b0001'1111) {
        // Long-form tag number
        tag.number = 0;
        while (offset < buffer.size()) {
            const Byte b = buffer[offset++];
            tag.number = (tag.number << 7) | (b & 0x7F);
            if ((b & 0x80) == 0) {
                break;
            }
        }
    }

    return tag;
}

struct BerLength {
    bool indefinite = false;
    std::size_t length = 0;
};

BerLength ReadLength(const std::vector<Byte> &buffer, std::size_t &offset)
{
    if (offset >= buffer.size()) {
        throw PinParseError("Unexpected end of buffer while reading BER length");
    }

    const Byte first = buffer[offset++];
    if (first == 0x80) {
        return {true, 0};
    }

    if ((first & 0x80) == 0) {
        return {false, first};
    }

    const std::uint8_t num_bytes = first & 0x7F;
    if (num_bytes == 0 || num_bytes > 8) {
        throw PinParseError("Unsupported BER length size");
    }

    std::size_t length = 0;
    for (std::uint8_t i = 0; i < num_bytes; ++i) {
        if (offset >= buffer.size()) {
            throw PinParseError("Unexpected end of buffer while reading BER length body");
        }
        length = (length << 8) | buffer[offset++];
    }

    return {false, length};
}

bool IsEoc(const std::vector<Byte> &buffer, std::size_t offset)
{
    return offset + 1 < buffer.size() && buffer[offset] == 0x00 && buffer[offset + 1] == 0x00;
}

void SkipElement(const std::vector<Byte> &buffer, std::size_t &offset)
{
    const BerTag tag = ReadTag(buffer, offset);
    const BerLength len = ReadLength(buffer, offset);

    if (len.indefinite) {
        if (!tag.constructed) {
            throw PinParseError("Indefinite length used with primitive element");
        }

        while (true) {
            if (IsEoc(buffer, offset)) {
                offset += 2;
                break;
            }
            SkipElement(buffer, offset);
        }
    } else {
        if (offset + len.length > buffer.size()) {
            throw PinParseError("BER element exceeds buffer size");
        }
        offset += len.length;
    }
}

std::int64_t ParseInteger(const std::vector<Byte> &buffer, std::size_t &offset, std::size_t length)
{
    if (length == 0 || offset + length > buffer.size()) {
        throw PinParseError("Invalid integer length");
    }

    std::int64_t value = (buffer[offset] & 0x80) ? -1 : 0;
    for (std::size_t i = 0; i < length; ++i) {
        value = (value << 8) | buffer[offset + i];
    }
    offset += length;
    return value;
}

std::string ParseString(const std::vector<Byte> &buffer, std::size_t &offset, std::size_t length)
{
    if (offset + length > buffer.size()) {
        throw PinParseError("String overruns buffer");
    }
    std::string result(reinterpret_cast<const char *>(&buffer[offset]), length);
    offset += length;
    return result;
}

std::uint64_t ReadLE64(const std::vector<Byte> &buffer, std::size_t &offset)
{
    if (offset + 8 > buffer.size()) {
        throw PinParseError("Unexpected end of file while reading 64-bit value");
    }

    std::uint64_t value = 0;
    for (int i = 7; i >= 0; --i) {
        value = (value << 8) | static_cast<std::uint64_t>(buffer[offset + i]);
    }
    offset += 8;
    return value;
}

std::string ReadPascalString(const std::vector<Byte> &buffer, std::size_t &offset)
{
    const std::uint32_t length = ReadBE32(buffer, offset);
    if (offset + length > buffer.size()) {
        throw PinParseError("String length exceeds file size");
    }

    std::string result(reinterpret_cast<const char *>(&buffer[offset]), length);
    offset += length;
    return result;
}

struct PinIndex {
    std::uint32_t version = 0;
    bool is_protein = false;
    std::uint32_t volume_number = 0; // only meaningful for version 5
    std::string title;
    std::string lmdb_file; // version 5 only
    std::string date;
    std::uint32_t num_oids = 0;
    std::uint64_t total_length = 0;
    std::uint32_t max_length = 0;
    std::vector<std::uint32_t> header_offsets;
    std::vector<std::uint32_t> sequence_offsets;
    std::vector<std::uint32_t> ambiguity_offsets; // nucleotide only
};

struct SeqId {
    std::string type;
    std::string value;
    std::optional<std::int64_t> version;
};

struct BlastDefLine {
    std::string title;
    std::vector<SeqId> seqids;
    std::optional<std::int64_t> taxid;
};

PinIndex ParsePinFile(const std::filesystem::path &path)
{
    const auto data = ReadFile(path);
    std::size_t offset = 0;

    PinIndex index;
    index.version = ReadBE32(data, offset);
    if (index.version != 4 && index.version != 5) {
        throw PinParseError("Unsupported database format version: " + std::to_string(index.version));
    }

    const std::uint32_t seq_type_flag = ReadBE32(data, offset);
    index.is_protein = (seq_type_flag == 1);

    if (index.version == 5) {
        index.volume_number = ReadBE32(data, offset);
    }

    index.title = ReadPascalString(data, offset);

    if (index.version == 5) {
        index.lmdb_file = ReadPascalString(data, offset);
    }

    index.date = ReadPascalString(data, offset);

    index.num_oids = ReadBE32(data, offset);
    index.total_length = ReadLE64(data, offset);
    index.max_length = ReadBE32(data, offset);

    const std::size_t count = static_cast<std::size_t>(index.num_oids) + 1;
    auto read_offset_array = [&](std::vector<std::uint32_t> &target) {
        target.resize(count);
        for (std::size_t i = 0; i < count; ++i) {
            target[i] = ReadBE32(data, offset);
        }
    };

    read_offset_array(index.header_offsets);
    read_offset_array(index.sequence_offsets);

    if (!index.is_protein) {
        read_offset_array(index.ambiguity_offsets);
    }

    if (offset != data.size()) {
        // The legacy format sometimes appends extra data; keep a small guard so we
        // at least warn that we ignored it.
        std::cerr << "Warning: trailing bytes in index file after parsing known fields." << std::endl;
    }

    return index;
}

std::filesystem::path DerivePhrPath(const std::filesystem::path &pin_path)
{
    auto result = pin_path;
    result.replace_extension(".phr");
    return result;
}

std::vector<std::string> ExtractHeaders(const PinIndex &index, const std::filesystem::path &phr_path)
{
    const auto data = ReadFile(phr_path);
    std::vector<std::string> headers;
    headers.reserve(index.num_oids);

    for (std::size_t i = 0; i < index.num_oids; ++i) {
        const std::uint32_t start = index.header_offsets.at(i);
        const std::uint32_t end = index.header_offsets.at(i + 1);
        if (end < start || end > data.size()) {
            throw std::runtime_error("Header offsets for OID " + std::to_string(i) + " are invalid");
        }
        headers.emplace_back(reinterpret_cast<const char *>(&data[start]), end - start);
    }

    return headers;
}

bool IsVisibleLikeTag(const BerTag &tag)
{
    // Blast titles are encoded as VisibleString in most databases, but older
    // volumes sometimes use different string types. Accept all universal
    // string encodings so we decode the title even when the concrete tag
    // varies (e.g., PrintableString instead of VisibleString).
    if (tag.cls != BerClass::Universal) {
        return false;
    }

    switch (tag.number) {
    case 12: // UTF8String
    case 18: // NumericString
    case 19: // PrintableString
    case 20: // TeletexString
    case 21: // VideotexString
    case 22: // IA5String
    case 25: // GraphicString
    case 26: // VisibleString
    case 27: // GeneralString
    case 28: // UniversalString
    case 29: // CharacterString
    case 30: // BMPString
        return true;
    default:
        return false;
    }
}

std::string ParseVisible(const std::vector<Byte> &buffer, std::size_t &offset)
{
    const BerTag inner_tag = ReadTag(buffer, offset);
    const BerLength inner_len = ReadLength(buffer, offset);
    if (!IsVisibleLikeTag(inner_tag)) {
        throw PinParseError("Expected string type inside explicit tag");
    }

    if (!inner_tag.constructed) {
        if (inner_len.indefinite) {
            throw PinParseError("Primitive string used with indefinite length");
        }
        return ParseString(buffer, offset, inner_len.length);
    }

    const bool inner_indef = inner_len.indefinite;
    const std::size_t inner_end = inner_indef ? buffer.size() : offset + inner_len.length;
    std::string combined;

    while (true) {
        if (inner_indef && IsEoc(buffer, offset)) {
            offset += 2;
            break;
        }
        if (!inner_indef && offset >= inner_end) {
            break;
        }

        const BerTag chunk_tag = ReadTag(buffer, offset);
        const BerLength chunk_len = ReadLength(buffer, offset);
        if (IsVisibleLikeTag(chunk_tag) && !chunk_tag.constructed && !chunk_len.indefinite) {
            combined += ParseString(buffer, offset, chunk_len.length);
        } else {
            SkipElement(buffer, offset);
        }
    }

    if (!inner_indef && offset < inner_end) {
        offset = inner_end;
    }

    return combined;
}

std::optional<std::string> ExtractVisibleLike(const std::vector<Byte> &buffer,
                                              std::size_t &offset, std::size_t limit)
{
    while (offset < limit) {
        if (IsEoc(buffer, offset)) {
            offset += 2;
            break;
        }

        const std::size_t element_start = offset;
        const BerTag tag = ReadTag(buffer, offset);
        const BerLength len = ReadLength(buffer, offset);

        if (IsVisibleLikeTag(tag)) {
            if (tag.constructed) {
                const bool indef = len.indefinite;
                const std::size_t end = indef ? limit : offset + len.length;
                if (auto inner = ExtractVisibleLike(buffer, offset, end)) {
                    return inner;
                }
                if (!indef && offset < end) {
                    offset = end;
                }
            } else if (len.indefinite) {
                throw PinParseError("Primitive string used with indefinite length");
            } else {
                return ParseString(buffer, offset, len.length);
            }
        } else if (len.indefinite) {
            if (!tag.constructed) {
                throw PinParseError("Indefinite length used with primitive element");
            }
            while (true) {
                if (IsEoc(buffer, offset)) {
                    offset += 2;
                    break;
                }
                if (auto inner = ExtractVisibleLike(buffer, offset, limit)) {
                    return inner;
                }
            }
        } else {
            if (offset + len.length > buffer.size()) {
                throw PinParseError("BER element exceeds buffer size");
            }
            offset += len.length;
        }

        if (offset <= element_start) {
            // Safety: ensure forward progress to avoid infinite loops if the
            // input is malformed.
            throw PinParseError("Failed to advance while scanning for string element");
        }
    }

    return std::nullopt;
}

std::string TagNameFromNumber(std::uint32_t num)
{
    static const std::unordered_map<std::uint32_t, std::string> kNames = {
        {0, "local"},        {1, "gibbsq"},      {2, "gibbmt"},    {3, "giim"},
        {4, "genbank"},      {5, "embl"},        {6, "pir"},       {7, "swissprot"},
        {8, "patent"},       {9, "other"},       {10, "general"},  {11, "gi"},
        {12, "ddbj"},        {13, "prf"},        {14, "pdb"},      {15, "tpg"},
        {16, "tpe"},        {17, "tpd"},        {18, "gpipe"},    {19, "named-annot-track"}
    };

    auto it = kNames.find(num);
    return it != kNames.end() ? it->second : ("unknown-" + std::to_string(num));
}

std::string ParseExplicitVisible(const std::vector<Byte> &buffer, std::size_t &offset, const BerLength &len);
std::int64_t ParseExplicitInteger(const std::vector<Byte> &buffer, std::size_t &offset, const BerLength &len);

SeqId ParseTextSeqId(const std::vector<Byte> &buffer, std::size_t &offset, std::optional<std::size_t> end_limit)
{
    SeqId id;
    const BerLength len = ReadLength(buffer, offset);
    const bool indefinite = len.indefinite;
    const std::size_t end = indefinite ? (end_limit ? *end_limit : buffer.size()) : offset + len.length;

    while (true) {
        if (indefinite && IsEoc(buffer, offset)) {
            offset += 2; // consume Textseq-id EOC
            break;
        }
        if (!indefinite && offset >= end) {
            break;
        }

        const BerTag tag = ReadTag(buffer, offset);
        const BerLength field_len = ReadLength(buffer, offset);

        switch (tag.number) {
        case 0:
            if (id.value.empty()) {
                if (tag.constructed || field_len.indefinite) {
                    id.value = ParseExplicitVisible(buffer, offset, field_len);
                } else {
                    id.value = ParseString(buffer, offset, field_len.length);
                }
            } else if (field_len.indefinite) {
                SkipElement(buffer, offset);
            } else {
                offset += field_len.length;
            }
            break;
        case 1:
            if (tag.constructed || field_len.indefinite) {
                id.value = ParseExplicitVisible(buffer, offset, field_len);
            } else {
                id.value = ParseString(buffer, offset, field_len.length);
            }
            break;
        case 3:
            if (tag.constructed || field_len.indefinite) {
                id.version = ParseExplicitInteger(buffer, offset, field_len);
            } else {
                id.version = ParseInteger(buffer, offset, field_len.length);
            }
            break;
        default:
            if (field_len.indefinite) {
                SkipElement(buffer, offset);
            } else {
                offset += field_len.length;
            }
            break;
        }
    }

    if (indefinite) {
        offset += 2; // consume EOC
    }

    return id;
}

SeqId ParseSeqId(const std::vector<Byte> &buffer, std::size_t &offset)
{
    const std::size_t seqid_start = offset;
    const BerTag tag = ReadTag(buffer, offset);
    SeqId id;
    id.type = TagNameFromNumber(tag.number);

    if (tag.cls != BerClass::ContextSpecific) {
        throw PinParseError("Seq-id uses unexpected tag class");
    }

    if (tag.constructed) {
        if (tag.number == 14) {
            // PDB-seq-id ::= SEQUENCE { mol VisibleString, chain INTEGER DEFAULT 32, ... }
            const BerLength len = ReadLength(buffer, offset);
            const bool indefinite = len.indefinite;
            const std::size_t end = indefinite ? buffer.size() : offset + len.length;
            while (indefinite ? !IsEoc(buffer, offset) : offset < end) {
                const BerTag field_tag = ReadTag(buffer, offset);
                const BerLength field_len = ReadLength(buffer, offset);
                if (field_len.indefinite) {
                    throw PinParseError("Indefinite length inside PDB-seq-id");
                }
                if (field_tag.cls == BerClass::Universal && field_tag.number == 26 && id.value.empty()) {
                    id.value = ParseString(buffer, offset, field_len.length);
                } else if (field_tag.cls == BerClass::Universal && field_tag.number == 2 && !id.version) {
                    id.version = ParseInteger(buffer, offset, field_len.length);
                } else {
                    offset += field_len.length;
                }
            }
            if (indefinite) {
                offset += 2;
            }
        } else {
            // Textseq-id, Giimport-id, Dbtag, etc.
            id = ParseTextSeqId(buffer, offset, std::nullopt);
            id.type = TagNameFromNumber(tag.number);
        }
    } else {
        // Primitive encodings (INTEGER based choices)
        const BerLength len = ReadLength(buffer, offset);
        if (len.indefinite) {
            throw PinParseError("Unexpected indefinite length for primitive Seq-id");
        }
        const std::int64_t val = ParseInteger(buffer, offset, len.length);
        id.value = std::to_string(val);
    }

    const std::size_t seqid_end = offset;
    if (id.value.empty()) {
        std::string best;
        std::string current;
        for (std::size_t i = seqid_start; i < seqid_end; ++i) {
            const char ch = static_cast<char>(buffer[i]);
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' || ch == '.') {
                current.push_back(ch);
            } else {
                if (current.size() > best.size()) {
                    best.swap(current);
                }
                current.clear();
            }
        }
        if (current.size() > best.size()) {
            best.swap(current);
        }
        if (!best.empty()) {
            id.value = best;
        }
    }

    return id;
}

std::vector<SeqId> ParseSeqIdList(const std::vector<Byte> &buffer, std::size_t &offset)
{
    const BerTag tag = ReadTag(buffer, offset);
    if (tag.cls != BerClass::Universal || tag.number != 16 || !tag.constructed) {
        throw PinParseError("Expected SEQUENCE for Seq-id list");
    }

    const BerLength len = ReadLength(buffer, offset);
    const bool indefinite = len.indefinite;
    const std::size_t end = indefinite ? buffer.size() : offset + len.length;

    std::vector<SeqId> ids;
    while (true) {
        if (indefinite && IsEoc(buffer, offset)) {
            offset += 2;
            break;
        }
        if (!indefinite && offset >= end) {
            break;
        }
        ids.push_back(ParseSeqId(buffer, offset));
    }
    return ids;
}

std::vector<SeqId> ParseSeqIdField(const std::vector<Byte> &buffer, std::size_t &offset)
{
    const BerLength len = ReadLength(buffer, offset);
    const std::size_t start = offset;

    auto ids = ParseSeqIdList(buffer, offset);

    if (len.indefinite) {
        while (!IsEoc(buffer, offset)) {
            SkipElement(buffer, offset);
        }
        offset += 2;
    } else {
        const std::size_t end = start + len.length;
        if (offset < end) {
            offset = end;
        }
    }

    return ids;
}

std::int64_t ParseExplicitInteger(const std::vector<Byte> &buffer, std::size_t &offset, const BerLength &len)
{
    const std::size_t start = offset;
    const BerTag inner_tag = ReadTag(buffer, offset);
    const BerLength inner_len = ReadLength(buffer, offset);
    if (inner_tag.cls != BerClass::Universal || inner_tag.number != 2 || inner_len.indefinite) {
        throw PinParseError("Expected INTEGER inside explicit wrapper");
    }
    const std::int64_t value = ParseInteger(buffer, offset, inner_len.length);

    if (len.indefinite) {
        while (!IsEoc(buffer, offset)) {
            SkipElement(buffer, offset);
        }
        offset += 2;
    } else {
        const std::size_t end = start + len.length;
        if (offset < end) {
            offset = end;
        }
    }

    return value;
}

std::string ParseExplicitVisible(const std::vector<Byte> &buffer, std::size_t &offset, const BerLength &len)
{
    const std::size_t start = offset;
    const std::size_t end = len.indefinite ? buffer.size() : start + len.length;
    std::string result;

    try {
        result = ParseVisible(buffer, offset);
    } catch (const PinParseError &) {
        // Fall back to a more permissive scan in case the explicit wrapper
        // contains additional layers or unexpected ordering before the
        // string we need.
        offset = start;
        auto recovered = ExtractVisibleLike(buffer, offset, end);
        result = recovered.value_or(std::string());
    }

    if (len.indefinite) {
        while (offset < end && !IsEoc(buffer, offset)) {
            SkipElement(buffer, offset);
        }
        if (IsEoc(buffer, offset)) {
            offset += 2;
        }
    } else if (offset < end) {
        offset = end; // skip any trailing explicit content we did not decode
    }

    return result;
}

std::vector<BlastDefLine> DecodeDeflineSet(const std::string &blob, std::string *error_out = nullptr)
{
    const std::vector<Byte> buffer(blob.begin(), blob.end());
    std::size_t offset = 0;

    const BerTag outer_tag = ReadTag(buffer, offset);
    if (outer_tag.cls != BerClass::Universal || outer_tag.number != 16) {
        throw PinParseError("Expected Blast-def-line-set sequence");
    }
    const BerLength outer_len = ReadLength(buffer, offset);
    const bool outer_indef = outer_len.indefinite;
    const std::size_t outer_end = outer_indef ? buffer.size() : offset + outer_len.length;

    std::vector<BlastDefLine> deflines;

    try {
        while (true) {
            if (outer_indef && IsEoc(buffer, offset)) {
                offset += 2;
                break;
            }
            if (!outer_indef && offset >= outer_end) {
                break;
            }

            const std::size_t def_start = offset;
            const BerTag def_tag = ReadTag(buffer, offset);
            if (def_tag.cls != BerClass::Universal || def_tag.number != 16 || !def_tag.constructed) {
                offset = def_start;
                SkipElement(buffer, offset);
                continue;
            }
            const BerLength def_len = ReadLength(buffer, offset);
            const bool def_indef = def_len.indefinite;
            const std::size_t def_end = def_indef ? buffer.size() : offset + def_len.length;

            BlastDefLine entry;

            try {
                while (true) {
                    if (def_indef && IsEoc(buffer, offset)) {
                        offset += 2;
                        break;
                    }
                    if (!def_indef && offset >= def_end) {
                        break;
                    }

                    const BerTag field_tag = ReadTag(buffer, offset);
                    if (field_tag.cls != BerClass::ContextSpecific) {
                        SkipElement(buffer, offset);
                        continue;
                    }

                    switch (field_tag.number) {
                    case 0: { // title
                        const BerLength len = ReadLength(buffer, offset);
                        if (field_tag.constructed || len.indefinite) {
                            entry.title = ParseExplicitVisible(buffer, offset, len);
                        } else {
                            entry.title = ParseString(buffer, offset, len.length);
                        }
                        break;
                    }
                    case 1: { // seqid list
                        entry.seqids = ParseSeqIdField(buffer, offset);
                        break;
                    }
                    case 2: { // taxid integer
                        const BerLength len = ReadLength(buffer, offset);
                        if (field_tag.constructed || len.indefinite) {
                            entry.taxid = ParseExplicitInteger(buffer, offset, len);
                        } else {
                            entry.taxid = ParseInteger(buffer, offset, len.length);
                        }
                        break;
                    }
                    default:
                        SkipElement(buffer, offset);
                        break;
                    }
                }

                deflines.push_back(entry);
            } catch (const std::exception &ex) {
                if (error_out && error_out->empty()) {
                    *error_out = ex.what();
                }
                if (!entry.title.empty() || !entry.seqids.empty() || entry.taxid) {
                    deflines.push_back(entry);
                }
                break;
            }
        }
    } catch (const std::exception &ex) {
        if (error_out) {
            *error_out = ex.what();
        }
    }

    return deflines;
}

void DumpHeaders(const std::vector<std::string> &headers, const std::filesystem::path &output_dir)
{
    std::filesystem::create_directories(output_dir);
    for (std::size_t i = 0; i < headers.size(); ++i) {
        const auto file_path = output_dir / ("header_" + std::to_string(i) + ".bin");
        std::ofstream out(file_path, std::ios::binary);
        out.write(headers[i].data(), static_cast<std::streamsize>(headers[i].size()));
    }
}

std::string TruncateForDisplay(const std::string &data, std::size_t max_bytes = 32)
{
    std::ostringstream os;
    os << std::hex << std::setfill('0');
    const std::size_t limit = std::min(max_bytes, data.size());
    for (std::size_t i = 0; i < limit; ++i) {
        os << std::setw(2) << (static_cast<unsigned>(static_cast<unsigned char>(data[i])));
        if (i + 1 != limit) {
            os << ' ';
        }
    }
    if (data.size() > limit) {
        os << " ...";
    }
    return os.str();
}

std::string FormatSeqId(const SeqId &id)
{
    std::ostringstream os;
    os << id.type << ':' << (id.value.empty() ? "<none>" : id.value);
    if (id.version) {
        os << '.' << *id.version;
    }
    return os.str();
}

void PrintSummary(const PinIndex &index, const std::vector<std::string> &headers)
{
    std::cout << "Database version : " << index.version << "\n";
    std::cout << "Sequence type    : " << (index.is_protein ? "protein" : "nucleotide") << "\n";
    std::cout << "Volume number    : " << index.volume_number << "\n";
    std::cout << "Title            : " << index.title << "\n";
    if (!index.lmdb_file.empty()) {
        std::cout << "LMDB file        : " << index.lmdb_file << "\n";
    }
    std::cout << "Date             : " << index.date << "\n";
    std::cout << "Sequences        : " << index.num_oids << "\n";
    std::cout << "Total length     : " << index.total_length << "\n";
    std::cout << "Max sequence len : " << index.max_length << "\n";
    std::cout << "\nHeader blocks:\n";

    for (std::size_t i = 0; i < headers.size(); ++i) {
        std::cout << "  OID " << i << " -> " << headers[i].size() << " bytes" << '\n';
        std::string decode_error;
        const auto deflines = DecodeDeflineSet(headers[i], &decode_error);
        if (deflines.empty()) {
            std::cout << "    (no deflines decoded)" << '\n';
        }
        for (std::size_t j = 0; j < deflines.size(); ++j) {
            const auto &def = deflines[j];
            std::cout << "    Defline " << j << ": "
                      << (def.title.empty() ? "<no title>" : def.title) << '\n';
            if (!def.seqids.empty()) {
                std::cout << "      IDs    : ";
                for (std::size_t k = 0; k < def.seqids.size(); ++k) {
                    if (k) {
                        std::cout << ", ";
                    }
                    std::cout << FormatSeqId(def.seqids[k]);
                }
                std::cout << '\n';
            }
            if (def.taxid) {
                std::cout << "      TaxID : " << *def.taxid << '\n';
            }
        }
        if (!decode_error.empty()) {
            std::cout << "    Warning: partial decode - " << decode_error << '\n';
            std::cout << "    Raw: " << TruncateForDisplay(headers[i]) << '\n';
        }
    }
}

} // namespace

int main(int argc, char *argv[])
{
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <database.pin> [output_dir]" << std::endl;
        return 1;
    }

    try {
        const std::filesystem::path pin_path = argv[1];
        const std::filesystem::path phr_path = DerivePhrPath(pin_path);

        const PinIndex index = ParsePinFile(pin_path);
        const auto headers = ExtractHeaders(index, phr_path);

        if (argc == 3) {
            DumpHeaders(headers, argv[2]);
        }

        PrintSummary(index, headers);
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 2;
    }

    return 0;
}
