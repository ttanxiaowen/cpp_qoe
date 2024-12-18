//
// Created by TXW on 24-12-15.
//

#ifndef RESOURCE_H
#define RESOURCE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <regex>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

using namespace std;

// 常量定义
constexpr double global_PI = 3.141592653589793;
constexpr double global_C_light = 299792457.999999984;  // 光速
constexpr double global_Re = 6378.1363e3;  // 地球半径（单位：米）
constexpr double global_K = -228.6;  // 玻尔兹曼常数
const unordered_set<string> facility_name = { "shanghai", "nanjing", "beijing", "hainan", "lasa", "ulu" };


struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

// CNR和FE表的定义
constexpr double CNR_table[] = {
    -2.35, -1.24, -0.3, 1.00, 2.23, 3.10, 4.03, 5.18, 6.2, 7.91, 9.35, 10.69,
    11.03, 12.89, 13.64, 14.28, 15.69, 16.05
};

constexpr double FE_table[] = {
    0.490243, 0.656448, 0.789412, 0.988858, 1.188304, 1.322253, 1.487473, 1.654663,
    1.766451, 2.22, 2.47, 2.64, 3.1, 3.5, 3.9, 4.119540, 4.3, 4.4
};
static double CNR2FrequenceEff(const double CNR) {
    if (CNR <= CNR_table[0])
        return 0;
    int  i;
    for (i = 0; i < 18; ++i)
    {
        if (CNR_table[i] > CNR)break;
    }
    return FE_table[i - 1]*1e6;
}


// 单个坐标的3D距离计算
static double Distance3D(const  vector<double>& a) {
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// 个坐标点之间的3D距离计算
static  double Distance3D(const  vector<double>& a, const  vector<double> &b) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) +
        (a[1] - b[1]) * (a[1] - b[1]) +
        (a[2] - b[2]) * (a[2] - b[2]));
}
// 辅助数据结构
class Data {
public:
    vector<double> p;
    vector<double> v;
    vector<double> lla;

    explicit Data(const vector<double>& p,
        const vector<double>& v = vector<double>(),
        const vector<double>& lla = vector<double>())
        : p(p), v(v), lla(lla) {
    }
    Data() = default;
};

class Traffic {
public:
    int id{};
    string src;
    string target;
    double rate{};

    // Constructor
    Traffic(const int id, string  src, string  target, const double rate)
        : id(id), src(std::move(src)), target(std::move(target)), rate(rate) {
    }

    Traffic() = default;
};

class Sat_link {
public:
    double total{};
    double remain{};

    explicit Sat_link(double total) : total(total), remain(total) {}

    void reset() {
        remain = total;
    }
    Sat_link() = default;
};

class Path {
public:
    vector<string> path;
    double rtt{};


    Path(const vector<string>& path, const double rtt) : path(path), rtt(rtt) {}
    Path() = default;
};

class SatAntenna {
public:
    double UsrLink_UL_G_t = 1.1;
    double UsrLink_UL_Frequency = 30 * 1e9;
    double UsrLink_UL_Bandwith = 50 * 1e6;
    double UsrLink_UL_BeamNum = 48;
    double UsrLink_DL_Frequency = 2000000000;
    double UsrLink_DL_Bandwith = 20 * 1e6;
    double UsrLink_DL_BeamNum = 48;
    double UserLink_EIRP_density = 34;  // dB W/MHz
    double UserLink_DL_Gt = 30;  // dBi

    double FeedLink_UL_G_t = 13;
    double FeedLink_UL_Frequency = 30 * 1e9;
    double FeedLink_UL_Bandwith = 100 * 1e6;
    double FeedLink_DL_Frequency = 20 * 1e9;
    double FeedLink_DL_Bandwith = 100 * 1e6;  // Hz
    double FeedLink_EIRP_density = 4;  // dB W/MHz
    double FeedLink_DL_Gt = 38.5;  // dBi

    double ISL_tx_capacity = 8000 * 1e6;
    double ISL_rx_capacity = 8000 * 1e6;
    string working_type = "TM";


};

class TerminalAntenna {
public:
    // 直接赋值
    double G_t = 34.2;  // 最大发送增益，单位：dB
    double G_r = 34.2;  // 最大接收增益，单位：dB
    double T = 150;     // 天线温度，单位：K
    double NF = 7;      // 噪声，单位：dB
    double f_UL = 20 * 1e9;  // 上行频点，单位���Hz
    double f_DL = 20 * 1e9;  // 下行频点，单位：Hz
    double Pt = 3.2;    // 发射机功率，单位：w
    double EIRP = 39.1; // 等效全向辐射功率，单位：dBW/MHz
    double G_T = 0;     // 接收系统性能因数，单位：dB/K
};

class FacilityAntenna {
public:
    // 直接赋值
    double G_t = 43.2;    // 最大发送增益，单位：dB
    double G_r = 39.7;    // 最大接收增益，单位：dB
    double BW_UL = 30 * 1e9; // 上行带宽分配，单位：Hz
    double BW_DL = 20 * 1e9; // 下行带宽分配，单位：Hz
    double T = 150;       // 天线温度，单位：K
    double NF = 1.2;      // 噪声，单位：dB
    double f_UL = 30 * 1e9; // 上行频点，单位：Hz
    double f_DL = 20 * 1e9; // 下行频点，单位：Hz
    double Pt = 2;        // 发射机功率，单位：w
    double EIRP = 0;      // 等效全向辐射功率，单位：dBW/MHz
    double G_T = 0;       // 接收系统性能因数，单位：dB/K
};
const SatAntenna sat_antenna;
const TerminalAntenna user_antenna;
const FacilityAntenna fac_antenna;
class Satellite {
public:
    int id{};
    string name;

    double total_user_up_capacity{};
    double total_user_down_capacity{};
    double total_faci_up_capacity{};
    double total_faci_down_capacity{};
    double total_use_power{};
    double total_faci_power{};

    double remaining_user_up_capacity{};
    double remaining_user_down_capacity{};
    double remaining_faci_up_capacity{};
    double remaining_faci_down_capacity{};
    double remaining_user_power{};
    double remaining_faci_power{};
    Satellite() = default;
    explicit Satellite(const vector<double>& p)
        :total_use_power(1000), total_faci_power(1000) {
        // 初始化
        this->total_user_up_capacity = cal_beam_center_up_link(p);
        this->total_user_down_capacity = cal_beam_center_down_link(p);
        this->total_faci_up_capacity = cal_faci_beam_center_up_link(p);
        this->total_faci_down_capacity = cal_faci_beam_center_down_link(p);
        this->remaining_user_up_capacity = total_user_up_capacity;
        this->remaining_user_down_capacity = total_user_down_capacity;
        this->remaining_faci_up_capacity = total_faci_up_capacity;
        this->remaining_faci_down_capacity = total_faci_down_capacity;
        this->remaining_user_power = total_use_power;
        this->remaining_faci_power = total_faci_power;
    }

    void reset() {
        remaining_user_up_capacity = total_user_up_capacity;
        remaining_user_down_capacity = total_user_down_capacity;
        remaining_faci_down_capacity = total_faci_down_capacity;
        remaining_faci_up_capacity = total_faci_up_capacity;
        remaining_user_power = total_use_power;
        remaining_faci_power = total_faci_power;
    }

    static double cal_beam_center_up_link(const vector<double>& sat) {


        // 计算自由空间损耗
        double Lp = 20 * log10(4 * global_PI * (Distance3D(sat) - global_Re) * sat_antenna.UsrLink_UL_Frequency / global_C_light);
        // 计算载噪比
        double C_N = user_antenna.EIRP + sat_antenna.UsrLink_UL_G_t - Lp - global_K - 10 * log10(sat_antenna.UsrLink_UL_Bandwith);
        return CNR2FrequenceEff(C_N) * sat_antenna.UsrLink_UL_Bandwith * 48/ 1e6;
    }

    static double cal_beam_center_down_link(const  vector<double>& sat) {
        // 计算EIRP
        double EIRP = pow(10, (sat_antenna.UserLink_EIRP_density / 10.0));
        EIRP = 10 * log10(EIRP * sat_antenna.UsrLink_DL_Bandwith / 1e6);
        // 计算自由空间损耗
        double Lp = 20 * log10(4 * global_PI * (Distance3D(sat) - global_Re) * sat_antenna.UsrLink_DL_Frequency / global_C_light);
        // 计算噪声功率
        double noise = 10 * log10(10) + 10 * log10(sat_antenna.UsrLink_DL_Bandwith)+ user_antenna.NF;
        // 计算载噪比
        double C_N = EIRP + user_antenna.G_r - Lp - global_K - noise;
        return CNR2FrequenceEff(C_N) * sat_antenna.UsrLink_DL_Bandwith * 48 / 1e6;
    }

    static double cal_faci_beam_center_down_link(const  vector<double>& sat) {
        // 计算EIRP
        double EIRP = pow(10, sat_antenna.FeedLink_EIRP_density / 10.0);
        EIRP = 10 * log10(EIRP * sat_antenna.FeedLink_DL_Bandwith / 1e6);
        // 计算自由空间损耗
        double Lp = 20 * log10(4 * global_PI * (Distance3D(sat) - global_Re) * sat_antenna.FeedLink_DL_Frequency / global_C_light);
        double noise = 10 * log10(fac_antenna.Pt) + fac_antenna.NF + 10 * log10(sat_antenna.FeedLink_DL_Bandwith);
        // 计算载噪比
        double C_N = EIRP + fac_antenna.G_r - Lp - global_K - noise;

        return CNR2FrequenceEff(C_N) * sat_antenna.FeedLink_DL_Bandwith *5/ 1e6;
    }

    // 计算上行链路的CNR和频率效率
    static double cal_faci_beam_center_up_link(const  vector<double>& sat) {
        // 计算自由空间损耗
        double Lp = 20 * log10(4 * global_PI * (Distance3D(sat) - global_Re) * sat_antenna.FeedLink_UL_Frequency / global_C_light);

        // 计算载噪比
        double C_N = 10 * log10(fac_antenna.Pt) + fac_antenna.G_t + sat_antenna.FeedLink_UL_G_t - Lp - global_K - 10 * log10(sat_antenna.FeedLink_UL_Bandwith);

        return CNR2FrequenceEff(C_N) * sat_antenna.FeedLink_UL_Bandwith *5/ 1e6;
    }


};

const unordered_map<string, Data> facility_position = {
    {"shanghai",Data({-2847751.2850821805, 4651738.5544014331,3306400.0246972283})},
    {"nanjing",Data({-2604132.68593206,4736841.10283461,3385497.98159804})},
    {"beijing",Data({-2176280.657815475, 4382058.474235897, 4091739.9482194087})},
    {"hainan",Data({-2036462.1147977512, 5672736.4802588429, 2087893.7795304032})},
    {"lasa",Data({-113522.47828731654, 5545270.592883083, 3156627.3137886818})},
    {"ulu",Data({191604.80967738523, 4598496.0708747767, 4416798.9852317292})},
};

class Info {
public:
    std::filesystem::path pathroot;
    std::vector<std::vector<std::vector<double>>> matrix;
    unordered_map<int, unordered_map<string, Data>> user_position;
    unordered_map<int, unordered_map<string, Data>> satellite_position;
    unordered_map<int, unordered_map<string, unordered_map<string, Sat_link>>> satellite_link;
    unordered_map<int, unordered_map<string, unordered_map<string, Sat_link>>> satellite_link_copy;
    unordered_map<string, int> sat_id_map;
    unordered_map<int, string> id_sat_map;
    unordered_map<int, unordered_map<string, Satellite>> satellite;
    unordered_map<int, unordered_map<string, Satellite>> satellite_copy;
    unordered_map<int, unordered_map<int, vector<Path>>> links;
    vector<Traffic> traffic;
    unordered_set<string> users;
    int num_users;
    int slot = -1;
    int num_time_slots;
    int num_actions = 20;
    int length = 0;

    void init(int input_length, int input_slot,std::filesystem::path start_path) {
        // 初始化
        pathroot =find_root_directory(start_path);
        num_time_slots = 31;
        num_users = input_length;
        slot = input_slot;
        length = input_length;  // 别忘了设置 length
        cout << "start" << endl;
        // 读取数据
        matrix = dijiskra_matrix_read();
        user_position = user_position_read();
        satellite_position = satellite_position_read();

        links = link_read();
        traffic_read(traffic, users);
        idmap_read(sat_id_map, id_sat_map, satellite);

        satellite_link = satellite_link_read();
        satellite_link_copy = satellite_link;
        satellite_copy = satellite;
        cout << "over" << endl;
    }

    static std::filesystem::path find_root_directory(std::filesystem::path start_path) {

        std::filesystem::path current_path = start_path;
        cout << "root" << current_path<< endl;
        while (current_path != current_path.root_path()) {
            // 假设项目根目录包含名为 CMakeLists.txt 的文件
            if (std::filesystem::exists(current_path / "CMakeLists.txt")) {
                return current_path;
            }
            current_path = current_path.parent_path();  // 返回上级目录
        }

        return current_path;  // 如果未找到根目录，返回当前目录
    }
    static void parse_vector_from_string(const string& str, const regex& regex, vector<double>& vec) {
        smatch match;
        if (regex_search(str, match, regex)) {
            string basic_string = match[1].str();
            stringstream ss(basic_string);
            string item;
            while (getline(ss, item, ',')) {
                erase_if(item, ::isspace);
                vec.push_back(stod(item)); // 将每个逗号分隔的子串转为 int
            }

        }
    }

    static vector<int> parse_points(const string& points_str) {
        vector<int> points;
        stringstream ss(points_str);
        string item;
        while (getline(ss, item, ',')) {
            points.push_back(stoi(item));
        }
        return points;
    }

    std::vector<std::vector<std::vector<double>>> dijiskra_matrix_read() {
        std::cout << "dijiskra_matrix_read……" << std::endl;

        std::vector<std::vector<std::vector<double>>> matrix;
        std::vector<std::vector<double>> current_layer;
        std::filesystem::path file_path = "data\\matrix.csv";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            std::cerr << "Failed to open file!" << std::endl;
            return matrix;
        }

        std::string line;
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t") + 1);

            if (line.empty()) {
                if (!current_layer.empty()) {
                    matrix.push_back(current_layer);
                    current_layer.clear();
                }
            }
            else {
                std::vector<std::string> row = split(line, ',');
                std::vector<double> current_row;
                for (const auto& item : row) {
                    try {
                        current_row.push_back(std::stof(item));
                    }
                    catch (const std::invalid_argument& e) {
                        std::cerr << "Invalid number format in row: " << item << ", error: " << e.what() << std::endl;
                    }
                }
                current_layer.push_back(current_row);
            }
        }
        if (!current_layer.empty()) {
            matrix.push_back(current_layer);
        }
        std::cout << "matrix_read已完成" << std::endl;
        return matrix;
    }


    unordered_map<int, unordered_map<string, Data>> user_position_read() {
        cout << "user_position_read……" << endl;
        unordered_map<int, unordered_map<string, Data>> user_position;
        std::filesystem::path file_path = "data\\terminal_p.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            cerr << "无法打开文件: terminal_p.txt" << endl;
            return user_position;
        }
        string line;
        regex time_regex(R"('Time':(\d+))");
        regex terminal_regex(R"('Terminal':(\w+))");
        regex lla_regex(R"('LLA':\[(.*?)\])"); // 提取 LLA
        regex p_regex(R"('P':\[(.*?)\])"); // 提取 P
        regex v_regex(R"('V':\[(.*?)\])"); // 提取 V

        while (getline(file, line)) {
            smatch match;
            if (regex_search(line, match, time_regex)) {
                int time = stoi(match[1]);
                if (time > num_time_slots) continue;
                if (slot != -1 && (time < slot || time > slot)) continue;

                if (regex_search(line, match, terminal_regex)) {
                    string Terminal = match[1];

                    vector<double> lla, p, v;
                    parse_vector_from_string(line, lla_regex, lla);

                    parse_vector_from_string(line, p_regex, p);
                    parse_vector_from_string(line, v_regex, v);


                    if (user_position.find(time) == user_position.end()) {
                        user_position[time] = {};
                    }
                    user_position[time][Terminal] = Data(p, v, lla);
                }
            }
        }
        cout << "user_position_read已完成" << endl;
        return user_position;
    }

    unordered_map<int, unordered_map<string, Data>> satellite_position_read() {
        cout << "satellite_position_read……" << endl;
        unordered_map<int, unordered_map<string, Data>> satellite_position;
        std::filesystem::path file_path = "data\\satellite_p.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            cerr << "无法打开文件: satellite_p.txt" << endl;
            return satellite_position;
        }
        string line;
        regex time_regex(R"('Time':(\d+))");
        regex sat_name_regex(R"('sat_name':(\w+))");
        regex p_regex(R"('P':\[(.*?)\])");
        regex v_regex(R"('V':\[(.*?)\])");

        while (getline(file, line)) {
            smatch match;
            if (regex_search(line, match, time_regex)) {
                int time = stoi(match[1]);
                if (time > num_time_slots) continue;
                if (slot != -1 && (time < slot || time > slot)) continue;

                if (regex_search(line, match, sat_name_regex)) {
                    string sat_name = match[1];

                    vector<double> p, v;
                    parse_vector_from_string(line, p_regex, p);
                    parse_vector_from_string(line, v_regex, v);

                    if (satellite_position.find(time) == satellite_position.end()) {
                        satellite_position[time] = {};
                    }
                    satellite_position[time][sat_name] = Data(p, v);
                }
            }
        }

        cout << "satellite_position_read已完成" << endl;
        return satellite_position;
    }

    unordered_map<int, unordered_map<string, unordered_map<string, Sat_link>>> satellite_link_read() {
        cout << "satellite_link_read……" << endl;
        unordered_map<int, unordered_map<string, unordered_map<string, Sat_link>>> satellite_link;

        std::filesystem::path file_path = "data\\satellite_link.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            cerr << "无法打开文件: satellite_link.txt" << endl;
            return satellite_link;
        }
        string line;
        SatAntenna sat_antenna;
        regex point_regex(R"(\[([^\]]+)\])");


        while (getline(file, line)) {
            smatch match;
            if (regex_search(line, match, regex(R"('Time':(\d+))"))) {
                int time = stoi(match[1]);
                if (time > num_time_slots) {
                    continue;
                }
                if (slot != -1 && (time < slot || time > slot)) {
                    continue;
                }

                if (regex_search(line, match, regex(R"('sat_name':(\w+))"))) {
                    string sat_name = match[1];
                    // 提取连接点
                    if (regex_search(line, match, point_regex)) {
                        string points_str = match[1].str();
                        vector<int> points = parse_points(points_str);

                        for (auto point : points) {
                            string sat = id_sat_map.at(point); // 使用at()来抛出异常如果键未找到
                            satellite_link[time][sat_name][sat] = Sat_link(sat_antenna.ISL_tx_capacity);
                        }
                    }
                }
            }
        }

        cout << "satellite_link_read已完成" << endl;
        return satellite_link;
    }

    void traffic_read(vector<Traffic>& traffic, unordered_set<string>& user) {
        cout << "traffic_read……" << endl;


        std::filesystem::path file_path = "data\\traffic.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            cerr << "无法打开文件: traffic.txt" << endl;
            return;
        }
        string line;
        int id = 0;

        // Define regex patterns
        regex src_regex(R"('src_name':(\w+))");
        regex target_regex(R"('target_name':(\w+))");
        regex rate_regex(R"('needrate':([\d.]+))");

        // Read the first 'length' lines from the file
        for (size_t i = 0; i < length && getline(file, line); ++i) {
            smatch match;

            // Extract src_name
            if (regex_search(line, match, src_regex)) {
                string src_name = match[1];

                // Extract target_name
                if (regex_search(line, match, target_regex)) {
                    string target_name = match[1];

                    // Extract needrate
                    if (regex_search(line, match, rate_regex)) {
                        double rate = stod(match[1]);

                        // Create Traffic object and add to vector
                        traffic.emplace_back(id, src_name, target_name, rate);
                        id++;

                        // Add to user set based on src_name condition
                        if (src_name.find("car") == string::npos) {
                            user.insert(target_name);
                        }
                        else {
                            user.insert(src_name);
                        }
                    }
                }
            }
        }

        cout << "traffic_read已完成" << endl;
    }

    void idmap_read(unordered_map<string, int>& sat_id_map, unordered_map<int, string>& id_sat_map, unordered_map<int, unordered_map<string, Satellite>>& satellite) {
        cout << "idmap_read……" << endl;

        std::filesystem::path file_path = "data\\satIdMap.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            cerr << "无法打开文件: satIdMap.txt" << endl;
            return;
        }

        string line;
        while (getline(file, line)) {
            std::istringstream iss(line);
            std::string key;
            int value;

            if (iss >> key >> value) {
                sat_id_map[key] = value;
                id_sat_map[value] = key;
            }


            if (facility_name.find(key) == facility_name.end() && value != -1) {
                if (slot != -1) {
                    int time = slot;
                    if (satellite.find(time) == satellite.end()) {
                        satellite[time] = {};
                    }
                    satellite[time][key] = Satellite(satellite_position.at(time).at(key).p);
                }
                else {
                    for (int time = 0; time < num_time_slots; ++time) {
                        if (satellite.find(time) == satellite.end()) {
                            satellite[time] = {};
                        }
                        satellite[time][key] = Satellite(satellite_position.at(time).at(key).p);
                    }
                }
            }
        }

        cout << "idmap_read已完成" << endl;
    }

    static vector<string> splitString( string& input, char delimiter) {
        vector<string> result;
        stringstream ss(input);
        string token;

        while (getline(ss, token, delimiter)) {
            result.push_back(token);
        }
        return result;
    }

    unordered_map<int, unordered_map<int, vector<Path>>> link_read()  {
        cout << "link_read……" << endl;

        std::filesystem::path file_path = "data\\upadate_link_1000.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);

        if (!file.is_open()) {
            cerr << "无法打开文件 upadate_link_1000.txt" << endl;
            return links;
        }
        unordered_map<int, unordered_map<int, vector<Path>>> links;
        string line;
        regex time_regex(R"('Time':(\d+))");
        regex traffic_id_regex(R"('traffic_id':(\w+))");
        regex path_regex(R"(\[([^\]]+)\])");
        regex rtt_regex(R"('rtt':([\d.]+))");

        while (getline(file, line)) {
            smatch match;

            // 提取 time
            if (regex_search(line, match, time_regex)) {
                int time = stoi(match[1].str());

                if (time > num_time_slots) continue;
                if (slot != -1 && time < slot) continue;
                if (slot != -1 && time > slot) continue;

                // 提取 traffic_id
                if (regex_search(line, match, traffic_id_regex)) {
                    int traffic_id = stoi(match[1].str());
                    if (traffic_id >= length) continue;

                    // 提取 path
                    if (regex_search(line, match, path_regex)) {
                        string path = match[1].str();
                        vector<string> satelliteVector = splitString(path, ',');
                        // 提取 rtt
                        if (regex_search(line, match, rtt_regex)) {
                            double rtt = stod(match[1].str());

                            // 确保字典结构
                            if (links.find(time) == links.end()) {
                                links[time] = {};
                            }
                            if (links[time].find(traffic_id) == links[time].end()) {
                                links[time][traffic_id] = {};
                            }

                            // 添加连接点
                            links[time][traffic_id].emplace_back(satelliteVector, rtt);
                        }
                    }
                }
            }
        }

        file.close();
        cout << "link_read已完成" << endl;
        return links;
    }

    unordered_map<int, unordered_map<string, vector<string>>> facility_link_read() {
        cout << "facility_link_read..." << endl;
        unordered_map<int, unordered_map<string, vector<string>>> facility_link;

        std::filesystem::path file_path = "data\\facility_link.txt";
        std::filesystem::path full_path = pathroot / file_path;
        ifstream file(full_path);
        if (!file.is_open()) {
            cerr << "无法打开文件: facility_link.txt." << endl;
            return facility_link;
        }

        string line;
        regex time_regex("'Time':(\\d+)");
        regex facility_regex("'facility':(\\w+)");
        regex points_regex("\\[([^]]+)\\]");

        while (getline(file, line)) {
            smatch match;

            // Parse 'Time'
            if (!regex_search(line, match, time_regex)) {
                continue;
            }
            int time = stoi(match[1].str());

            if (time > num_time_slots) {
                continue;
            }
            if (slot != -1 && time < slot) {
                continue;
            }
            if (slot != -1 && time > slot) {
                continue;
            }

            // Parse 'facility'
            if (!regex_search(line, match, facility_regex)) {
                continue;
            }
            string facility = match[1].str();

            // Parse connection points
            if (!regex_search(line, match, points_regex)) {
                continue;
            }
            string points_str = match[1].str();
            vector<string> connection_points = split(points_str, ',');

            // Ensure dictionary structure
            facility_link[time][facility] = connection_points;
        }

        cout << "facility_link_read complete." << endl;
        return facility_link;
    }

    static vector<string> split(const string& str, char delimiter) {
        vector<string> tokens;
        stringstream ss(str);
        string token;
        while (getline(ss, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }
};


#endif //RESOURCE_H
