#include <algorithm>
#include <boost/spirit/include/qi.hpp>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

using namespace BamTools;

namespace Cigar {
    enum class OpType : char {
        MATCH       = 'M',
        INS         = 'I',
        DEL         = 'D',
        SKIP        = 'N',
        SOFTCLIP    = 'S',
        HARDCLIP    = 'H',
        PADDING     = 'P',
        SEQMATCH    = '=',
        SEQMISMATCH = 'X',
        AMBIGUOUS   = 'A'
    };

    struct Op {
        OpType type;
        unsigned len;
    };

    typedef std::vector<Op> String;

    std::ostream& operator<<(std::ostream& s, Op const& x);
    std::ostream& operator<<(std::ostream& s, String const& x);
    bool operator!=(Op const& x, Op const& y);

    template<typename Container> String fromMdTag(Container const& x);
    template<typename Container> String fromCigarOps(Container const& x);

    // Implementation

    template<typename Iterator>
    class MdTokenizer_ {
    public:
        typedef typename std::iterator_traits<Iterator>::value_type value_type;

        MdTokenizer_(Iterator beg, Iterator end) : beg_(beg), end_(end) {}

        bool eof() const {
            return beg_ == end_;
        }

        static bool test(value_type x) {
            return !std::isalpha(x);
        }

        // Don't call this when this->eof() is true
        Op nextOperation() {
            assert(!eof());
            namespace qi = boost::spirit::qi;

            Op op{OpType::SEQMISMATCH, 0u};

            // qi::parse advances beg_ when it successfully parses
            if (qi::parse(beg_, end_, qi::uint_, op.len)) {
                op.type = OpType::SEQMATCH;
            } else {
                if (*beg_ == '^') {
                    ++beg_; // Don't count caret towards length
                    op.type = OpType::DEL;
                }

                auto where = std::find_if(beg_, end_, test);

                assert(std::distance(beg_, where) > 0);
                op.len = unsigned(std::distance(beg_, where));
                beg_ = where;
            }

            return op;
        }

    private:
        Iterator beg_;
        Iterator end_;
    };

    template<typename Container>
    String fromMdTag(Container const& x) {
        auto beg = x.begin();
        auto end = x.end();
        MdTokenizer_<decltype(beg)> tokenizer(beg, end);
        String result;
        while (!tokenizer.eof()) {
            Op next = tokenizer.nextOperation();
            if (!(next.type == OpType::SEQMATCH && next.len == 0)) {
                result.push_back(next);
            }
        }
        return result;
    }

    template<typename Container>
    String fromCigarOps(Container const& x) {
        String result;
        for (auto op = x.begin(); op != x.end(); op++) {
            result.push_back(Op {static_cast<OpType>(op->Type), op->Length});
        }
        return result;
    }

    std::ostream& operator<<(std::ostream& s, Op const& x) {
        s << x.len << static_cast<char>(x.type);
        return s;
    }

    std::ostream& operator<<(std::ostream& s, String const& x) {
        std::copy(x.begin(), x.end(), std::ostream_iterator<Op>(s));
        return s;
    }

    bool operator!=(Op const& x, Op const& y) {
        return (x.len != y.len) || (x.type != y.type);
    }
}

class Result {
public:
    Result(BamAlignment & al) : _alignment(al) {
        if(!_alignment.BuildCharData())
            throw std::runtime_error("Error unpacking alignment strings");
    }
    // TODO do we need other things like copy, move constructor or destructors?

    Cigar::String compose_ops() {
        try {
            return _compose_ops();
            // TODO we would return here
        } catch (std::exception& e) {
            std::string name;
            if(_alignment.BuildCharData()) {
                name = "[" + _alignment.Name + "] ";
            } else {
                name = "<unknown> ";
            }

            throw std::runtime_error("\n[31m" + name + e.what() + "[0m");
        }
    }

private:
    BamAlignment & _alignment;
    std::string _query_str;
    Cigar::String _cigar_str;
    Cigar::String _md_str;
    Cigar::String _comp_str;
    Cigar::String::iterator _cigar;
    Cigar::String::iterator _cigar_end;
    Cigar::String::iterator _md;
    Cigar::String::iterator _md_end;
    std::string::iterator _query;
    std::string::iterator _query_end;

    Cigar::String _compose_ops() {
        _cigar_str = Cigar::fromCigarOps(_alignment.CigarData);
        _md_str    = Cigar::fromMdTag(_extract_md());
        _comp_str  = Cigar::String();
        _query_str = _alignment.QueryBases;

        _cigar     = _cigar_str.begin();
        _cigar_end = _cigar_str.end();
        _md        = _md_str.begin();
        _md_end    = _md_str.end();
        _query     = _query_str.begin();
        _query_end = _query_str.end();

        while(_cigar != _cigar_end) {
            iterate();
        }

        if (_md != _md_end)
            throw std::runtime_error("Did not reach end of MD tag");
        if (_query != _query_end)
            throw std::runtime_error("Did not reach end of query sequence");

        return _comp_str; // TODO should copy out, right?
    }

    void iterate() {
        switch (_cigar->type) {
            using Cigar::OpType;

            case OpType::HARDCLIP:
                hardclip();
                break;
            case OpType::SOFTCLIP:
                softclip();
                break;
            case OpType::INS:
                insertion();
                break;
            case OpType::DEL:
                deletion();
                break;
            case OpType::MATCH:
            case OpType::SEQMATCH:
            case OpType::SEQMISMATCH:
                match();
                break;
            case OpType::SKIP:
                skip();
                break;
            case OpType::PADDING:
                padding();
                break;
            default:
                throw std::runtime_error("Invalid Cigar string");
        }
    }

    void hardclip() {
        record(*_cigar);
        _cigar++;
    }

    void softclip() {
        record(*_cigar);
        _query += _cigar->len;
        _cigar++;
    }

    void insertion() {
        record(*_cigar);
        _query += _cigar->len;
        _cigar++;
    }

    void deletion() {
        if (_md >= _md_end)
            throw std::runtime_error("Tried to dereference MD iterator past end of MD tag");
        if (*_md != *_cigar)
            throw std::runtime_error("Cigar indicated deletion but MD did not match");

        record(*_cigar);

        _cigar++;
        _md++;
    }

    void match() {
        if (_md >= _md_end)
            throw std::runtime_error("Tried to dereference MD iterator past end of MD tag");
        if (_query >= _query_end)
            throw std::runtime_error("Tried to dereference query iterator past end of query sequence");

        auto & md = *_md;
        auto & cigar = *_cigar;
        auto & query = *_query;

        if ((cigar.type != Cigar::OpType::MATCH) && (cigar.type != md.type))
            throw std::runtime_error("Cigar and MD don't have compatible (mis)match types");

        if ( query == 'N' || query == 'n' ) {
            cigar.len--;
            md.len--;

            _query++;

            if (cigar.len == 0) {
                _cigar++;
            }
            if (md.len == 0) {
                _md++;
            }

            record(Cigar::Op {Cigar::OpType::AMBIGUOUS, 1u});

            return;
        }

        if (cigar.len < md.len) {
            record(Cigar::Op {md.type, cigar.len});
            md.len -= cigar.len;
            _cigar++;
            _query += cigar.len;
        } else if (cigar.len > md.len) {
            record(md);
            cigar.len -= md.len;
            _md++;
            _query += md.len;
        } else {
            record(md);
            _cigar++;
            _md++;
            _query += md.len;
        }
    }

    void skip() {
        throw std::runtime_error("Skip (S) in Cigar unimplemented");
    }

    void padding() {
        throw std::runtime_error("Padding (P) in Cigar unimplemented");
    }

    inline void record(Cigar::Op const& op) {
        // TODO How is memory working here? Is there only one copy?
        _comp_str.push_back(op);
    }

    std::string _extract_md() {
        std::string destination;
        if (!_alignment.HasTag("MD"))
            throw std::runtime_error("Expected MD tag");
        if (!_alignment.GetTag("MD", destination))
            throw std::runtime_error("Error getting MD tag");
        return destination;
    }
};

namespace Count {
    template<typename item_type>
    class Single {
    public:
        void add_count(uint64_t start, uint64_t length) {
            uint64_t end = start + length;
            while (end > counts.size()) {
                counts.push_back(0u);
            }

            for (uint64_t index = start; index < end; index++) {
                counts[index] += 1;
            }
        }

        item_type operator[](uint64_t index) {
            return index >= counts.size() ? 0u : counts[index];
        }

        uint64_t size() {
            return counts.size();
        }

        item_type sum() {
            item_type sum = 0;
            for (auto val = counts.begin(); val != counts.end(); val++) {
                sum += *val;
            }
            return sum;
        }

    private:
        std::vector<item_type> counts;
    };

    typedef Count::Single<uint64_t> Counter;

    class All {
    public:
        All(uint32_t read_end) : read_end(read_end), _used(false) {}

        void process_cigar(Cigar::String cigar, bool reverse) {
            _used = true;
            if (reverse) {
                _process_cigar(cigar.rbegin(), cigar.rend());
            } else {
                _process_cigar(cigar.begin(), cigar.end());
            }
        }

        uint64_t max_size() {
            uint64_t max = 0;

            max = std::max(insertion.size(), max);
            max = std::max(deletion.size(), max);
            max = std::max(ambiguous.size(), max);
            max = std::max(mismatch.size(), max);
            max = std::max(match.size(), max);
            max = std::max(total.size(), max);

            return max;
        }

        bool used() {
            return _used;
        }

        std::string table() {
            uint64_t max = max_size();
            std::stringstream table;

            for (uint64_t i = 0; i < max; i++) {
                table << position(i);
            }

            table << sums();

            return table.str();
        }

        const char * header() {
            return
                "read_end\t"
                "position\t"
                "total\t"
                "match\t"
                "error\t"
                "error_rate\t"
                "mismatch\t"
                "mismatch_rate\t"
                "ambiguous\t"
                "ambiguous_rate\t"
                "insertion\t"
                "insertion_rate\t"
                "deletion\t"
                "deletion_rate\n";
        }

    private:
        Counter insertion;
        Counter deletion;
        Counter ambiguous;
        Counter mismatch;
        Counter match;
        Counter total;
        uint32_t read_end;
        bool _used;

        template<typename iterator_type>
        void _process_cigar(iterator_type start, iterator_type end) {
            uint64_t position = 0u;

            while (start != end) {
                position = _process_op(*start, position);
                start++;
            }
        }

        uint64_t _process_op(Cigar::Op const & op, uint64_t position) {
            switch (op.type) {
                using Cigar::OpType;

                case OpType::HARDCLIP:
                    return position + op.len;
                case OpType::SOFTCLIP:
                    return position + op.len;
                case OpType::INS:
                    insertion.add_count( position, op.len );
                    total.add_count( position, op.len );
                    return position + op.len;
                case OpType::DEL:
                    for (unsigned times = 0; times < op.len; times++) {
                        deletion.add_count( position - 1, 1 ); // TODO verify this
                    }
                    return position;
                case OpType::SEQMATCH:
                    match.add_count( position, op.len );
                    total.add_count( position, op.len );
                    return position + op.len;
                case OpType::SEQMISMATCH:
                    mismatch.add_count( position, op.len );
                    total.add_count( position, op.len );
                    return position + op.len;
                case OpType::AMBIGUOUS:
                    ambiguous.add_count( position, op.len );
                    total.add_count( position, op.len );
                    return position + op.len;
                default:
                    throw std::runtime_error("Invalid CigarOp for Counts");
            }
        }

        std::string position(uint64_t position) {
            std::stringstream line;

            auto total_     = total[position];
            auto match_     = match[position];
            auto mismatch_  = mismatch[position];
            auto ambiguous_ = ambiguous[position];
            auto insertion_ = insertion[position];
            auto deletion_  = deletion[position];

            line << read_end << "\t"
                << position + 1 << "\t"
                << tallies(total_, match_, mismatch_, ambiguous_, insertion_, deletion_)
                << std::endl;

            return line.str();
        }

        std::string sums() {
            std::stringstream line;

            auto total_     = total.sum();
            auto match_     = match.sum();
            auto mismatch_  = mismatch.sum();
            auto ambiguous_ = ambiguous.sum();
            auto insertion_ = insertion.sum();
            auto deletion_  = deletion.sum();

            line << read_end << "\t"
                << "SUM" << "\t"
                << tallies(total_, match_, mismatch_, ambiguous_, insertion_, deletion_)
                << std::endl;

            return line.str();
        }

        std::string tallies(
                uint64_t total_, uint64_t match_,
                uint64_t mismatch_, uint64_t ambiguous_,
                uint64_t insertion_, uint64_t deletion_) {

            auto error_     = mismatch_ + ambiguous_ + insertion_ + deletion_;

            std::stringstream tallies;
            tallies.setf(std::ios::fixed, std::ios::floatfield);
            tallies.precision(3);

            tallies << total_ << "\t"
                << match_ << "\t"
                << error_ << "\t"
                << (double) error_ / total_ << "\t"
                << mismatch_ << "\t"
                << (double) mismatch_ / total_ << "\t"
                << ambiguous_ << "\t"
                << (double) ambiguous_ / total_ << "\t"
                << insertion_ << "\t"
                << (double) insertion_ / total_ << "\t"
                << deletion_ << "\t"
                << (double) deletion_ / total_;

            return tallies.str();
        }
    };
}


void process_alignment(BamAlignment &al, Count::All & single, Count::All & first, Count::All & second) {
    uint32_t flag = al.AlignmentFlag;
    if ( (flag & 0x4) || (flag & 0x100) || (flag & 0x400) ) {
        return;
    }

    Result results(al);
    Cigar::String composed = results.compose_ops();

    bool reverse = al.IsReverseStrand();

    if (al.IsFirstMate()) {
        first.process_cigar(composed, reverse);
    } else if (al.IsSecondMate()) {
        second.process_cigar(composed, reverse);
    } else {
        single.process_cigar(composed, reverse);
    }
}

uint64_t process_file(std::string const& path, Count::All & single, Count::All & first, Count::All & second) {
    uint64_t count = 0;

    BamReader reader;

    if (!reader.Open(path)) {
        throw std::invalid_argument("unable to open " + path);
    }

    BamAlignment al;
    while (reader.GetNextAlignmentCore(al)) {
        process_alignment(al, single, first, second);
        count++;
    }

    reader.Close();

    return count;
}

int main(int argc, char* argv[]) {
    uint64_t total = 0;

    if (argc < 2) {
        std::cerr << "Please provide at least one bam\n";
        return 1;
    }

    std::vector<std::string> inputFilenames;
    for (int i = 1; i < argc; i++) {
        inputFilenames.push_back(std::string(argv[i]));
    }

    Count::All single(0);
    Count::All first(1);
    Count::All second(2);

    for (auto input = inputFilenames.begin(); input != inputFilenames.end(); input++) {
        uint64_t count = process_file(*input, single, first, second);
        // std::cerr << count << "\t" << input << std::endl;
        total += count;
    }
    // std::cerr << total << "\t*" << std::endl;
    std::cout << single.header(); // Only print the header once
    if (single.used()) {
        std::cout << single.table();
    }
    if (first.used()) {
        std::cout << first.table();
    }
    if (second.used()) {
        std::cout << second.table();
    }
}


