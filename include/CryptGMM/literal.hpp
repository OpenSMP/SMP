#ifndef CRYPTCONV_LITERAL_HPP
#define CRYPTCONV_LITERAL_HPP
#include <vector>
#include <string>

std::string trim(const std::string &line);
std::vector<std::string> splitBySpace(const std::string &line);
std::pair<double, double> mean_std(std::vector<double> const& times);

#endif // CRYPTCONV_LITERAL_HPP
