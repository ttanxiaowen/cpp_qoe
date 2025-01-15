// //
// // Created by TTTxw on 2025/1/7.
// //

// #include "min.h"
// void caculate_min(int length,filesystem::path current_path) {
//     Info info;
//     info.init(length, -1,current_path);
//     const int ROUND = 6;
//     string suanfa = "MINDELAY";
//     string path = "max" ;
//     string file = path + "\\"+ to_string(length) + "result.txt";


//     float QoE = 0.0;

//     ofstream fileout("min" + to_string(length) + ".txt");
//     fileout << "num 算法, Time, traffic, rate, delay, loss, Qoe, select_path\n";

//     for (int time_slot = 0; time_slot < info.num_time_slots; ++time_slot) {
//         Dijskra dijskra;
//         auto cur_matrix = info.matrix[time_slot];
//         dijskra.initGraph(cur_matrix);

//         // 重置卫星链接
//         for (auto& sat : info.satellite_link[time_slot]) {
//             for (auto& link : sat.second) {
//                 link.second.reset();
//             }
//         }

//         // 处理每个流量
//         for (int i = 0; i < info.traffic.size(); ++i) {
//             auto src = info.traffic[i].src;
//             auto target = info.traffic[i].target;
//             bool f = false;

//             if (facility_name.find(src) == facility_name.end()) {
//                 vector<pair<int, float>> path_distance;

//                 if (info.user_canlink[time_slot][src].empty()) {
//                     cout << "traffic:" << i << " 没有可接星" << endl;
//                     // traffic_fail[i]++;
//                     continue;
//                 }

//                 for (auto sat : info.user_canlink[time_slot][src]) {
//                     double distance = Distance3D(info.user_position[time_slot][src].p, info.satellite_position[time_slot][sat].p);
//                     path_distance.push_back({sat, distance});
//                 }

//                 sort(path_distance.begin(), path_distance.end(), [](const pair<int, float>& a, const pair<int, float>& b) {
//                     return a.second < b.second;
//                 });

//                 for (const auto& [sat, dis] : path_distance) {
//                     double rtt = info.calculate_min_link()

//                     float rtt = dijskra.dijskra(info.sat_id_map[src], info.sat_id_map[target]);
//                     rtt += dis;
//                     vector<int> satellite_ids;
//                     dijskra.findWaypoint(satellite_ids, info.sat_id_map[src], info.sat_id_map[target]);
//                     vector<string> selected_link;
//                     for (const auto& id : satellite_ids) {
//                         selected_link.insert(info.id_sat_map(satellite_ids))
//                     }

//                     double rate = 0, loss = 0;
//                     vector<int> fail = {-2, -2};
//                     bool flag = can_meet_link_requirements(time_slot, i, selected_link, info, &rate, &loss, fail, -1);

//                     if (flag) {
//                         update_link(time_slot, selected_link, rate, info);
//                         float reward = calculate_qoe(rate, info.traffic[i].rate, rtt * 2 / global_C_light, loss);

//                         selected_link.insert(v.begin(), src);

//                         fileout << length << "\t" << suanfa << "\t" << time_slot << "\t" << info.traffic[i].id
//                                 << "\t" << need_capacity / info.traffic[i].rate << "\t" << rtt * 2 / global_C_light
//                                 << "\t" << loss << "\t" << reward << "\t" << join(satellite_ids, ",") << endl;

//                         f = true;
//                         break;
//                     }
//                 }

//                 if (!f) {
//                     traffic_fail[i]++;
//                     traffic_loss[i]++;
//                     fileout << length << " " << suanfa << " " << time_slot << " " << info.traffic[i].id
//                             << " 0 inf 1 0" << endl;
//                     cout << "没有接入星" << endl;
//                 }
//             } else {
//                 // 地面到用户的处理
//                 vector<pair<int, float>> path_distance;
//                 cout << "src:" << src << ", time:" << time_slot << endl;

//                 if (info.user_canlink[time_slot][target].empty()) {
//                     traffic_fail[i]++;
//                     continue;
//                 }

//                 for (int sat : info.user_canlink[time_slot][target]) {
//                     float distance = Distance3D(info.user_position[time_slot][target].p, info.satellite_position[time_slot][sat].p).distance();
//                     path_distance.push_back({sat, distance});
//                 }

//                 sort(path_distance.begin(), path_distance.end(), [](const pair<int, float>& a, const pair<int, float>& b) {
//                     return a.second < b.second;
//                 });

//                 for (const auto& [sat, dis] : path_distance) {
//                     float rtt = dijskra.dijskra(info.sat_id_map[src], info.sat_id_map[sat]);
//                     rtt += dis;
//                     vector<int> satellite_ids;
//                     dijskra.findWaypoint(satellite_ids, info.sat_id_map[src], info.sat_id_map[sat]);
//                     satellite_ids.push_back(target);

//                     auto [flag, need_capacity, loss, fail, unsatis, snr] = can_meet_link_requirements(time_slot, i, satellite_ids, info);

//                     if (flag) {
//                         update_link(time_slot, satellite_ids, need_capacity, info);
//                         float reward = calculate_qoe(need_capacity, info.traffic[i].rate, rtt * 2 / global_C_light, loss);
//                         user_qoe[i] += reward;
//                         userrtt[i] += rtt * 2 / global_C_light;
//                         rate_capacity[i] += need_capacity;
//                         QoE += reward;
//                         traffic_loss[i] += loss;
//                         user_snr.push_back(snr);

//                         fileout << length << "\t" << suanfa << "\t" << time_slot << "\t" << info.traffic[i].id
//                                 << "\t" << need_capacity / info.traffic[i].rate << "\t" << rtt * 2 / global_C_light
//                                 << "\t" << loss << "\t" << reward << "\t" << join(satellite_ids, ",") << endl;

//                         f = true;
//                         break;
//                     }
//                 }

//                 if (!f) {
//                     traffic_fail[i]++;
//                     traffic_loss[i]++;
//                     fileout << length << " " << suanfa << " " << time_slot << " " << info.traffic[i].id
//                             << " 0 inf 1 0" << endl;
//                 }
//             }
//         }
//     }

//     fileout << "Total QoE: " << QoE << endl;

//     // 处理输出
//     ofstream file("min_" + to_string(length) + ".txt");
//     file << "num, 算法, traffic, fail, Qoe, rtt, rate_ratio, loss\n";
//     for (int traffic_id = 0; traffic_id < info.traffic.size(); ++traffic_id) {
//         if (traffic_fail[traffic_id] != info.num_time_slots) {
//             file << length << " " << suanfa << " " << traffic_id << " " << traffic_fail[traffic_id] << " "
//                  << user_qoe[traffic_id] / info.num_time_slots << " "
//                  << userrtt[traffic_id] / (info.num_time_slots - traffic_fail[traffic_id]) << " "
//                  << rate_capacity[traffic_id] / info.num_time_slots / info.traffic[traffic_id].rate << " "
//                  << traffic_loss[traffic_id] / info.num_time_slots << endl;
//         } else {
//             file << length << " " << suanfa << " " << traffic_id << " " << traffic_fail[traffic_id] << " "
//                  << "0 inf 0 1" << endl;
//         }
//     }

//     ofstream snrfile("snr.txt");
//     for (float snr : user_snr) {
//         snrfile << snr << endl;
//     }
// }

// int main(int argc, char* argv[]) {
//     if (argc > 1) {
//         try {
//             int number = stoi(argv[1]);
//             caculate_min(number, 0);
//         } catch (const exception& e) {
//             caculate_min(1000, 0);
//         }
//     } else {
//         cout << "No number provided. Please run the script with a number as an argument." << endl;
//         caculate_min(100, 0);
//     }

//     return 0;
// }