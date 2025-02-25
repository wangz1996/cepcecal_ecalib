#ifndef LOGGER_HH
#define LOGGER_HH
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <unordered_map>
#include <vector>

template<typename T>
struct TypeConverter {
    static std::string toString(const T& value) {
        std::ostringstream oss;
        oss << value;
        return oss.str(); 
    }
};

class Logger{
public:
    static Logger& getInstance(){
        static Logger instance;
        return instance;
    }

    template<typename T>
    void log(const std::string& key, const T& value){
        umap_log[key].emplace_back(TypeConverter<T>::toString(value));
    }

    void init(const std::string& fname){
        file_.open(fname, std::ios::out | std::ios::trunc);
    }
private:
    Logger() = default;
    ~Logger(){
        for(auto& [key, value] : umap_log){
            file_ << key;
            for(auto& val : value){
                file_ << " " << val;
            }
            file_ << std::endl;
        }
        file_.close();
    }
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    std::unordered_map<std::string,std::vector<std::string>> umap_log;
    std::ofstream file_;

};
#endif