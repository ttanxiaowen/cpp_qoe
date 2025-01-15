#include "max.h"

// 环境类，模拟环境中的状态
class Environment_max
{
public:
    Info info;
    int current_time_slot = 0;
    vector<vector<double>> cur_matrix;

    Environment_max(Info info) : info(info)
    {
    }

    void reset(int time_slot)
    {
        current_time_slot = time_slot;
        cur_matrix = info.matrix[current_time_slot];
        info.satellite[current_time_slot] = info.satellite_copy[current_time_slot];
        info.satellite_link[current_time_slot] = info.satellite_link_copy[current_time_slot];
    }

    pair<int, int> get_state(int traffic_id, int len)
    {
        return make_pair(traffic_id, len); // 返回状态元组
    }

    vector<Path> get_available_links(int traffic_id)
    {
        return info.links[current_time_slot][traffic_id]; // 返回可用的链路
    }

    double get_reward(double rate, double needRate, double delay, double loss)
    {
        return calculate_qoe(rate, needRate, delay, loss); // 获取奖励
    }

    vector<vector<double>> update_fail_matrix(const vector<vector<double>> cur_matrix, int i, int j)
    {
        // 创建 cur_matrix 的副本
        vector<vector<double>> matrix(cur_matrix);
        matrix[i][j] = INFINITY;
        matrix[j][i] = INFINITY; // 如果矩阵是对称的，也要更新对角元素

        return matrix;
    }

    double calculate_min_link(const vector<vector<double>> &cur_matrix, const string &src, const string &target,
                              vector<string> &path)
    {
        path.clear();
        Dijskra dijkstra;
        dijkstra.initGraph(cur_matrix);
        int s = info.sat_id_map[src];
        int t = info.sat_id_map[target];
        double distance = dijkstra.dijskra(s, t);
        if (distance == INFINITY)
        {
            return -1;
        }
        vector<int> path_id;
        dijkstra.findWaypoint(path_id, s, t);
        for (int node : path_id)
        {
            path.push_back(info.id_sat_map[node]);
        }
        return distance;
    }

    bool add_user_link(int cur_id, vector<int> &fail, vector<string> &satellite_ids)
    {
        bool flag = false;
        int round = 8;
        auto up_matrix = cur_matrix;
        string src = info.traffic[cur_id].src;
        string tar = info.traffic[cur_id].target;
        double rtt = 0;
        for (int step = 0; step < round && !flag; ++step)
        {
            int i = -1, j = -1;
            if (fail[0] == -1)
            {
                i = info.sat_id_map[src];
                j = info.sat_id_map[satellite_ids[fail[1]]];
            }
            else
            {
                i = info.sat_id_map[satellite_ids[fail[0]]];
                j = info.sat_id_map[satellite_ids[fail[1]]];
            }
            up_matrix = update_fail_matrix(up_matrix, i, j);
            if (src.find("car") != string::npos)
            {
                // Car 做traffic:
                rtt = calculate_min_link(up_matrix, satellite_ids[0], satellite_ids.back(), satellite_ids);
                if (rtt == -1)
                    break;
                rtt += Distance3D(info.satellite_position[current_time_slot][satellite_ids[0]].p,
                                  info.user_position[current_time_slot][src].p);
            }
            else
            {
                // 信关站 做traffic:
                rtt = calculate_min_link(up_matrix, src, satellite_ids[satellite_ids.size() - 2], satellite_ids);
                if (rtt == -1)
                    break;
                rtt += Distance3D(info.satellite_position[current_time_slot][satellite_ids[satellite_ids.size() - 1]].p,
                                  info.user_position[current_time_slot][tar].p);
                satellite_ids.erase(satellite_ids.begin());
                satellite_ids.push_back(tar);
            }
            double rate = 0, loss = 0;
            flag = can_meet_link_requirements2(current_time_slot, cur_id, satellite_ids, info, &rate, &loss, fail, -1);
        }

        if (flag)
        {
            cur_matrix = up_matrix;
            // Append the path to the links map
            info.links[current_time_slot][info.traffic[cur_id].id].emplace_back(satellite_ids, rtt);
            cout << "Traffic:" << info.traffic[cur_id].id << " 添加一条链路" << endl;
            return true;
        }
        else
        {
            int i = -1, j = -1;
            if (fail[0] == -1)
            {
                i = info.sat_id_map[src];
                j = info.sat_id_map[satellite_ids[fail[1]]];
            }
            else
            {
                i = info.sat_id_map[satellite_ids[fail[0]]];
                j = info.sat_id_map[satellite_ids[fail[1]]];
            }
            up_matrix = update_fail_matrix(up_matrix, i, j);
            cur_matrix = up_matrix;
            cout << "Traffic:" << info.traffic[cur_id].id << " 添加链路失败" << endl;
            return false;
        }

        
    }
};

// Q-learning 代理类
class QLearningAgent_max
{
public:
    int traffic_length;
    int num_actions;
    double learning_rate;
    double discount_factor;
    double exploration_rate;
    unordered_map<pair<int, int>, vector<double>, pair_hash> q_table; // Q表
    // vector<int
    std::vector<std::vector<std::vector<double>>> q_table1;

    unordered_map<pair<int, int>, int, pair_hash> max_values; // 最终结果
    unordered_map<pair<int, int>, int, pair_hash> fail_nums;

    explicit QLearningAgent_max(int traffic_length, int num_actions)
        : traffic_length(traffic_length), num_actions(num_actions), learning_rate(0.1), discount_factor(0.9),
          exploration_rate(1.0)
    {
        q_table1.resize(traffic_length); // Adjust based on the number of possible traffic_ids
        for (int i = 0; i < traffic_length; ++i)
        {
            q_table1[i].resize(num_actions); // Adjust based on the number of possible len values
            for (int j = 0; j < num_actions; ++j)
            {
                q_table1[i][j].resize(num_actions, 0.0); // Initialize actions for each state to 0
            }
        }
    }

    int choose_action(pair<int, int> state, int len)
    {
        if (rand() / double(RAND_MAX) < exploration_rate && exploration_rate > 0.02)
        {
            return rand() % min(len, num_actions); // 探索模式
        }
        else
        {
            auto it = q_table.find(state);
            if (it != q_table.end())
            {
                return distance(it->second.begin(), max_element(it->second.begin(), it->second.end())); // 利用模式
            }
            return 0; // 如果没有找到状态，随机选择
        }
    }

    void update(pair<int, int> state, int action, double reward, pair<int, int> next_state)
    {
        if (q_table.find(state) == q_table.end())
        {
            q_table[state] = vector<double>(num_actions, 0.0); // 初始化 Q 表
        }

        auto &q_values = q_table[state];
        double best_next_q = 0.0;
        if (q_table.find(next_state) != q_table.end())
        {
            best_next_q = *max_element(q_table[next_state].begin(), q_table[next_state].end()); // 获取下一个状态的最大 Q 值
        }

        double td_target = reward + discount_factor * best_next_q;
        q_values[action] += learning_rate * (td_target - q_values[action]); // 更新 Q 值
    }

    void calculate_max_q()
    {
        for (const auto &pair : q_table)
        {
            const auto &key = pair.first;
            const auto &values = pair.second;
            if (!values.empty())
            {
                double max_value = *max_element(values.begin(), values.end());
                size_t max_index = distance(values.begin(), find(values.begin(), values.end(), max_value));
                max_values[key] = max_index;
            }
        }
    }
};

// 处理每个time槽的函数
void process_time_slot_max(int length, int time_slot, Info info, const int episodes, double lossrate)
{
    // 排序
    auto compareByRate = [](const Traffic &a, const Traffic &b)
    {
        return a.rate < b.rate;
    };
    std::sort(info.traffic.begin(), info.traffic.end(), compareByRate);

    // 初始化
    Environment_max env(info);
    QLearningAgent_max agent(info.num_actions, info.num_actions);

    for (int episode = 0; episode < episodes; ++episode)
    {
        env.reset(time_slot);
        cout << "length" << length << "time" << time_slot << "episode:" << episode << endl;
        for (int i = 0; i < info.traffic.size(); ++i)
        {
            auto traffic = info.traffic[i];
            string src = traffic.src;
            string target = traffic.target;
            int traffic_id = traffic.id;
            double need_rate = traffic.rate;
            int len = env.get_available_links(traffic_id).size();
            if (len == 0)
                continue;
            auto state = env.get_state(traffic_id, len);
            int action = agent.choose_action(state, len);
            double reward = -5;
            // 获取链路，计算奖励
            auto available_links = env.get_available_links(traffic_id);
            // if (traffic_id == 213 ) {
            //     action = 0;
            //  }
            if (action < len)
            {
                auto link = available_links[action];
                auto selected_link = link.path;
                auto rtt = available_links[action].rtt;
                double rate, loss;
                vector<int> fail = {-2, -2};
                bool flag = can_meet_link_requirements2(time_slot, i, selected_link, env.info, &rate, &loss, fail,
                                                        episode);
                if (flag)
                {
                    update_link(time_slot, selected_link, rate, env.info);
                    // reward = env.get_reward(rate, need_rate, rtt * 2 / global_C_light, loss);
                    reward = rate / 1e3;
                }
                else
                {
                    if (len < env.info.num_actions && agent.fail_nums[state] < 700)
                    {
                        if ((fail[1] == 0 && info.users.find(src) != info.users.end()) || (fail[1] == selected_link.size() - 1 && info.users.find(src) != info.users.end()))
                        {
                        }
                        else
                        {
                            if (!env.add_user_link(i, fail, selected_link))
                            {
                                agent.fail_nums[state]++;
                            }
                        }
                    }
                }
            }
            auto next_state = env.get_state(i + 1, len);
            agent.update(state, action, reward, next_state);
        }

        agent.exploration_rate = max(0.01, agent.exploration_rate * lossrate); // 减少探索率
    }

    // 计算最大值
    agent.calculate_max_q();
    
    // 最终结果
    string path = "MAX2";
    if (!filesystem::exists(path))
    {
        filesystem::create_directory(path);
    }

    string file = path + "/" + to_string(length) + "_" + to_string(time_slot) + ".txt";

    ofstream outFile(file);

    string allfile = path + "/allmax.txt";

    ofstream alloutFile(allfile, ios::app);
    env.reset(time_slot);
    double allULRate = 0;
    double allDLRate = 0;
    double allLoss = 0;
    double allDelay = 0;
    double success = 0;
    double allRatio = 0;
    double allqoe = 0;

    for (int i = 0; i < length; ++i)
    {
        auto traffic = info.traffic[i];
        string src = traffic.src;
        string target = traffic.target;
        int traffic_id = traffic.id;
        double need_rate = traffic.rate;

        int len = env.get_available_links(traffic_id).size();
        auto state = env.get_state(traffic_id, len);
        int action = agent.max_values[state];
        double reward = 0;
        ;
        // 获取链路，计算奖励
        auto available_links = env.get_available_links(traffic_id);
        if (len > 0 && action < len)
        {
            auto link = available_links[action];
            auto selected_link = link.path;
            auto rtt = available_links[action].rtt;
            double rate, loss;
            vector<int> fail = {-2, -2};
            bool flag = can_meet_link_requirements2(time_slot, i, selected_link, env.info, &rate, &loss, fail,
                                                    episodes);
            if (flag)
            {
                double delay = rtt * 2 / global_C_light;
                update_link(time_slot, selected_link, rate, env.info);
                reward = env.get_reward(rate, need_rate, delay, loss);

                selected_link.insert(selected_link.begin(), src);
                outFile << length << "\t" << "MAX\t" << time_slot << "\t" << traffic_id << "\t" << rate / need_rate << "\t"
                        << delay << "\t" << loss << "\t" << reward << "\t" << joinWithCommas(selected_link) << endl;
                allRatio += rate / need_rate;
                allDelay += delay;
                allLoss += loss;
                allULRate += facility_name.contains(src) ? 0 : rate;
                allDLRate += facility_name.contains(src) ? rate : 0;
                success++;
                allqoe += reward;
            }
            else
            {
                outFile << length << "\t" << "MAX " << time_slot << "\t" << traffic_id << "\t" << 0 << "\t"
                        << "inf" << "\t" << "inf" << "\t" << 0 << "\t" << endl;
            }
        }
        else
        {
            outFile << length << "\t" << "MAX\t" << time_slot << "\t" << traffic_id << "\t" << 0 << "\t"
                    << "inf" << "\t" << "inf" << "\t" << 0 << "\t" << endl;
        }
    }
    // qoe loss delay ratio fail  ul dl
    alloutFile << length << "\t" << "MAX\t" << time_slot << "\t" << allqoe / length << "\t" << allLoss / success << "\t" << allDelay / success << "\t" << allRatio / success << "\t" << (length - success) / length << "\t" << allULRate << "\t" << allDLRate << endl;
    cout << allqoe << endl;
    outFile.close();
    alloutFile.close();
}

void process_time_slot_threaded_max(int length, int time_slot, Info info, const int episodes, double lossrate)
{
    process_time_slot_max(length, time_slot, info, episodes, lossrate);
}

void appendFileContent_max(const filesystem::path &sourceFile, std::ofstream &targetFile)
{
    std::ifstream inputFile(sourceFile);
    if (!inputFile.is_open())
    {
        std::cerr << "无法打开文件: " << sourceFile << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line))
    {
        targetFile << line << std::endl; // 将每一行内容写入目标文件
    }

    inputFile.close();
}

void calculate_all_max(filesystem::path current_path, int len,int choice,int plus)
{
    for (int length = len; length < len + plus && length <= 2000; length += 50)
    {
        int episodes = (length >= 1500) ? 13000 : (length >= 800) ? 10000
                                              : (length > 600)    ? 6000
                                                                  : 4000;
        double lossrate = (length >= 800) ? 0.9996 : 0.9995;
        Info info;
        info.init(length, -1, current_path); // 假设info->init函数已经定义
        double times[4][2] = {{0, 10}, {10, 20}, {20, 31},{0,31}};
        double range[2][2] = {{0,3},{3,4}};
        for (int i = range[choice][0]; i <range[choice][1]; i++)
        {
            auto start = chrono::high_resolution_clock::now();
            vector<thread> threads;
            for (int time = times[i][0]; time < times[i][1]; ++time)
                threads.emplace_back(process_time_slot_max, length, time, info, episodes, lossrate);

            for (auto &t : threads)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = end - start;

            cout << "Time measured: " << elapsed.count() << " seconds.\n";
        }
    }
}