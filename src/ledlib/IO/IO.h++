#include <string>

namespace LED::IO {

class Output {

public:
    /* If filename given, write to file; for empty filename write to screen */
    void InitOutput(const std::string& filename, const std::string& headline);

    /* Add line to standard output with three numbers */
    void AddToOutput(const double n1, const double n2, const double n3);

    /* Add line to standard output with two numbers */
    void AddToOutput2(const double n1, const double n2);

private:
    bool create_directory_if_needed(const std::string& filepath) const;
    std::string fOutputFilePath;
};

} // namespace LED::IO