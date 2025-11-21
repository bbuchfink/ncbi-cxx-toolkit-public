#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct PinIndex {
    int version = 0;
    int sequence_type = 0;  // 1 for protein in BLAST databases
    uint32_t volume_number = 0;
    std::string title;
    std::string lmdb_name;
    std::string creation_date;
    uint32_t num_sequences = 0;
    uint64_t total_residues = 0;
    uint32_t max_length = 0;
    std::vector<uint32_t> sequence_offsets;
};

class BinaryCursor {
public:
    explicit BinaryCursor(std::vector<uint8_t> data) : m_Data(std::move(data)) {}

    uint32_t ReadUInt32() {
        EnsureAvailable(4);
        uint32_t value = (static_cast<uint32_t>(m_Data[m_Pos]) << 24) |
                         (static_cast<uint32_t>(m_Data[m_Pos + 1]) << 16) |
                         (static_cast<uint32_t>(m_Data[m_Pos + 2]) << 8) |
                         (static_cast<uint32_t>(m_Data[m_Pos + 3]));
        m_Pos += 4;
        return value;
    }

    uint64_t ReadUInt64() {
        EnsureAvailable(8);
        uint64_t value = 0;
        for (int i = 0; i < 8; ++i) {
            value = (value << 8) | m_Data[m_Pos + i];
        }
        m_Pos += 8;
        return value;
    }

    std::string ReadString() {
        const uint32_t length = ReadUInt32();
        EnsureAvailable(length);
        std::string value(reinterpret_cast<const char*>(m_Data.data() + m_Pos),
                          length);
        m_Pos += length;
        return value;
    }

    size_t Position() const { return m_Pos; }

    const std::vector<uint8_t>& Data() const { return m_Data; }

private:
    void EnsureAvailable(size_t bytes) const {
        if (m_Pos + bytes > m_Data.size()) {
            throw std::runtime_error("PIN file ended unexpectedly");
        }
    }

    std::vector<uint8_t> m_Data;
    size_t m_Pos = 0;
};

std::vector<uint8_t> ReadFile(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Unable to open file: " + path);
    }
    std::vector<uint8_t> data((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    if (data.empty()) {
        throw std::runtime_error("File is empty: " + path);
    }
    return data;
}

PinIndex ParsePin(const std::string& path) {
    BinaryCursor cursor(ReadFile(path));

    PinIndex index;
    index.version = static_cast<int>(cursor.ReadUInt32());
    index.sequence_type = static_cast<int>(cursor.ReadUInt32());

    if (index.version != 4 && index.version != 5) {
        throw std::runtime_error("Unsupported PIN format version: " +
                                 std::to_string(index.version));
    }

    if (index.version == 5) {
        index.volume_number = cursor.ReadUInt32();
    }

    index.title = cursor.ReadString();
    if (index.version == 5) {
        index.lmdb_name = cursor.ReadString();
    }
    index.creation_date = cursor.ReadString();

    index.num_sequences = cursor.ReadUInt32();
    index.total_residues = cursor.ReadUInt64();
    index.max_length = cursor.ReadUInt32();

    const size_t region_bytes = static_cast<size_t>(index.num_sequences + 1) * 4;
    const size_t header_offsets = cursor.Position();
    const size_t sequence_offsets = header_offsets + region_bytes;

    if (sequence_offsets + region_bytes > cursor.Data().size()) {
        throw std::runtime_error("PIN offset tables are incomplete");
    }

    index.sequence_offsets.reserve(index.num_sequences + 1);
    for (uint32_t i = 0; i <= index.num_sequences; ++i) {
        const size_t pos = sequence_offsets + i * 4;
        const uint32_t value = (static_cast<uint32_t>(cursor.Data()[pos]) << 24) |
                               (static_cast<uint32_t>(cursor.Data()[pos + 1]) << 16) |
                               (static_cast<uint32_t>(cursor.Data()[pos + 2]) << 8) |
                               (static_cast<uint32_t>(cursor.Data()[pos + 3]));
        index.sequence_offsets.push_back(value);
    }

    if (index.sequence_offsets.size() < 2 ||
        index.sequence_offsets.front() >= index.sequence_offsets.back()) {
        throw std::runtime_error("PIN sequence offsets appear to be corrupt");
    }

    if (index.sequence_type != 1) {
        throw std::runtime_error("This reader only supports protein databases (type 1)");
    }

    return index;
}

char DecodeResidue(uint8_t code) {
    static const char kNcbistdaaToAscii[] = {
        '*', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
        'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'Z', 'U', 'O',
        'J', '-'};

    if (code == 0) {
        return '\0';
    }

    if (code < sizeof(kNcbistdaaToAscii)) {
        return kNcbistdaaToAscii[code];
    }

    return '?';
}

std::string DecodeSequence(const std::vector<uint8_t>& data, uint32_t start, uint32_t end) {
    if (start > end || end > data.size()) {
        throw std::runtime_error("Sequence offsets exceed PSQ file length");
    }

    std::string decoded;
    decoded.reserve(end - start);
    for (uint32_t pos = start; pos < end; ++pos) {
        const char aa = DecodeResidue(data[pos]);
        if (aa == '\0') {
            break;  // sequences are NUL-terminated in protein volumes
        }
        decoded.push_back(aa);
    }
    return decoded;
}

std::string StripExtension(const std::string& path, const std::string& ext) {
    if (path.size() >= ext.size() &&
        path.compare(path.size() - ext.size(), ext.size(), ext) == 0) {
        return path.substr(0, path.size() - ext.size());
    }
    return path;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: psq_reader <database path without extension or .pin/.psq path>\n";
        return 1;
    }

    std::string input = argv[1];
    input = StripExtension(StripExtension(input, ".pin"), ".psq");
    const std::string pin_path = input + ".pin";
    const std::string psq_path = input + ".psq";

    try {
        const PinIndex index = ParsePin(pin_path);
        const std::vector<uint8_t> psq_bytes = ReadFile(psq_path);

        std::cout << "PIN metadata\n";
        std::cout << "  Version: " << index.version << "\n";
        std::cout << "  Database type: protein\n";
        std::cout << "  Title: " << index.title << "\n";
        if (!index.lmdb_name.empty()) {
            std::cout << "  LMDB backing file: " << index.lmdb_name << "\n";
        }
        std::cout << "  Created: " << index.creation_date << "\n";
        std::cout << "  Sequences: " << index.num_sequences << "\n";
        std::cout << "  Total residues: " << index.total_residues << "\n";
        std::cout << "  Longest sequence: " << index.max_length << " residues\n";

        for (uint32_t i = 0; i < index.num_sequences; ++i) {
            const uint32_t start = index.sequence_offsets[i];
            const uint32_t end = index.sequence_offsets[i + 1];
            const std::string sequence = DecodeSequence(psq_bytes, start, end);
            std::cout << ">oid_" << i << " length=" << sequence.size() << "\n";
            std::cout << sequence << "\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
