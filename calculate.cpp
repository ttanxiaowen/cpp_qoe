//
// Created by TXW on 24-12-15.
//

#include "calculate.h"
#include "resource.h"

double play_time_comparison(double x0, double x1) {
    x0 = min(x0, 1.0);
    double coefficient_x0 = -176.0910098652557;
    double coefficient_x1 = 0.040544440750675745;
    double coefficient_x0_squared = 137.76092720848308;
    double coefficient_x0_x1 = -0.036297203251921074;
    double coefficient_x1_squared = -6.54659592466234e-06;
    double coefficient_x0_cubed = -35.037032776120924;
    double coefficient_x0_squared_x1 = 0.010292844426660324;
    double coefficient_x0_x1_squared = -1.2985954444286255e-06;
    double coefficient_x1_cubed = 4.468507494581608e-09;
    double intercept = 74.66995667558035;

    return (
        intercept
        + coefficient_x0 * x0
        + coefficient_x1 * x1
        + coefficient_x0_squared * x0 * x0
        + coefficient_x0_x1 * x0 * x1
        + coefficient_x1_squared * x1 * x1
        + coefficient_x0_cubed * x0 * x0 * x0
        + coefficient_x0_squared_x1 * x0 * x0 * x1
        + coefficient_x0_x1_squared * x0 * x1 * x1
        + coefficient_x1_cubed * x1 * x1 * x1
        );
}

double exp_loss(double x) {
    return 1.01 * exp(-5.12 * x);
}

double exp_time(double x) {
    return 5.1 * exp(-0.05 * (x - 1));
}

double custom_function(double x) {
    double m = 1 / (1 + exp(-0.3 * (x - 1)));
    return m;
}

double calculate_qoe(double rate, double need_rate, double delay, double loss) {
    // Calculate QoE_loss
    double QoE_loss = 5 * exp_loss(loss * 100);
    // Calculate QoE_time
    double QoE_time = play_time_comparison(rate / need_rate, delay * 1000);
    QoE_time = exp_time(QoE_time);  // Ensure custom_function is defined
    // Calculate overall QoE
    double weight_loss_adjusted = (QoE_loss < QoE_time) ? 0.75 : 0.5;
    double weight_time_adjusted = 1 - weight_loss_adjusted;
    double QoE = weight_loss_adjusted * QoE_loss + weight_time_adjusted * QoE_time;

    return QoE / 5.036678606218123 * 5;
}

bool can_meet_link_requirements(int time, int index, const vector<string>& selected_link, Info& info, double* rate, double* loss, vector<int>& unsatis) {
    double need_rate = info.traffic[index].rate;
    const  string src = info.traffic[index].src;
    string target = info.traffic[index].target;
    double snr = 0;

    for (int i = 0; i < selected_link.size(); ++i) {
        if (i == 0) {
            // 用户-卫星 上行
            if (facility_name.find(src) == facility_name.end()) {
                // 计算用户-卫星距离
                double distance = Distance3D(info.user_position[time][src].p, info.satellite_position[time][selected_link[i]].p);
                snr = calculate_user_up_SNR(distance);
                if (snr < 0) {
                    cout << "user:" << src << "与卫星" << selected_link[i] << "的snr" << snr << " 小于0 \n";
                    unsatis[0] = -1;
                    unsatis[1] = 0;
                    return false;
                }
                *loss = calculate_loss(snr);
                double calcu_rate = calculate_rate(snr, sat_antenna.UsrLink_UL_Bandwith);
                if (calcu_rate > need_rate) {
                    *rate = need_rate;
                }
                else {
                    *rate = calcu_rate;
                    unsatis[0] = -1;
                    unsatis[1] = 0;
                }
                // 判断链路是否满足要求
                if (info.satellite[time][selected_link[i]].remaining_user_up_capacity < *rate) {
                    if (info.satellite[time][selected_link[i]].remaining_user_up_capacity < *rate / 3) {
                        cout << "Time:" << time << ", Traffic" << index << "," << selected_link[i] << "用户链路上行容量不足\n";
                        unsatis[0] = -1;
                        unsatis[1] = 0;
                        return false;
                    }
                    else {
                        *rate = info.satellite[time][selected_link[i]].remaining_user_up_capacity;
                    }
                }
            }
            else {
                // 地面站-卫星 上行
                double distance = Distance3D(facility_position.at(src).p, info.satellite_position[time][selected_link[i]].p);
                snr = calculate_user_up_SNR(distance); // 假设用相同的计算
                double calcu_rate = calculate_rate(snr, sat_antenna.FeedLink_UL_Bandwith);
                if (calcu_rate > need_rate) {
                    *rate = need_rate;
                }
                else {
                    *rate = calcu_rate;
                    unsatis[0] = -1;
                    unsatis[1] = 0;
                }

                // 判断链路是否满足要求
                if (info.satellite[time][selected_link[i]].remaining_faci_up_capacity < *rate) {
                    cout << "Time:" << time << ", Traffic" << index << "," << selected_link[i] << "馈电链路上行容量不足\n";
                    unsatis[0] = -1;
                    unsatis[1] = 0;
                    return false;
                }
            }
        }
        else {
            if (info.users.contains(selected_link[i])) {
                double distance = Distance3D(info.user_position[time][selected_link[i]].p,info.satellite_position[time].at(selected_link[i - 1]).p);
                double user_snr = calculate_user_up_SNR(distance);
                double user_loss = calculate_loss(user_snr);
                double need_power = calculate_sat_user_need_power(*rate, sat_antenna.UsrLink_DL_Bandwith, distance);
                *loss = max(*loss, user_loss);

                if (info.satellite[time].at(selected_link[i - 1]).remaining_user_power < need_power) {
                    cout << "Traffic" << index << ", " << selected_link[i - 1]
                        << " 用户链路功率 " << info.satellite[time].at(selected_link[i - 1]).remaining_user_power
                        << " < " << need_power << "\n";
                    unsatis[0] = i - 1;
                    unsatis[1] = i;
                    return false;
                }
                if (info.satellite[time].at(selected_link[i - 1]).remaining_user_down_capacity < *rate) {
                    cout << "Traffic" << index << ", " << selected_link[i - 1]
                        << " 用户链路容量 " << info.satellite[time].at(selected_link[i - 1]).remaining_user_down_capacity
                        << " < " << need_rate << "\n";
                    unsatis[0] = i - 1;
                    unsatis[1] = i;
                    return false;
                }
            }
            // 星间链路
            else if (!facility_name.contains(selected_link[i]) && !facility_name.contains(selected_link[i - 1])) {
                // 如果是星间链路
                if (info.satellite_link[time].at(selected_link[i - 1])[selected_link[i]].remain < *rate) {
                    cout << "Time:" << time << ", Traffic" << index << "," << selected_link[i - 1] << "," << selected_link[i]
                        << "星间链路容量不足\n";
                    unsatis[0] = i - 1;
                    unsatis[1] = i;
                    return false;
                }
            }
            // 馈电下行
            else if (!facility_name.contains(selected_link[i - 1]) &&
                     facility_name.contains(selected_link[i])) {
                if (info.satellite[time].at(selected_link[i - 1]).remaining_faci_down_capacity < *rate) {
                    cout << "Time:" << time << ", Traffic" << index << "," << selected_link[i - 1] << "馈电链路下行容量不足\n";
                    unsatis[0] = i - 1;
                    unsatis[1] = i;
                    return false;
                }

                // 计算功率需求
                double distance = Distance3D(facility_position.at(selected_link[i]).p, info.satellite_position[time].at(selected_link[i - 1]).p);
                double need_power = calculate_sat_faci_need_power(*rate, sat_antenna.FeedLink_UL_Bandwith,distance);
                if (info.satellite[time].at(selected_link[i - 1]).remaining_faci_power < need_power) {
                    cout << "Time:" << time << ", Traffic" << index << "," << selected_link[i - 1] << "馈电链路功率不足\n";
                    unsatis[0] = i - 1;
                    unsatis[1] = i;
                    return false;
                }
            }
            // 馈电上行
            else if (!facility_name.contains(selected_link[i]) &&
                     facility_name.contains(selected_link[i - 1])) {

                if (info.satellite[time][selected_link[i]].remaining_faci_up_capacity < *rate) {
                    cout << "Time:" << time << ", Traffic" << index << "," << selected_link[i] << "馈电链路上行容量不足\n";
                    unsatis[0] = i - 1;
                    unsatis[1] = i;
                    return false;
                }
            }
            // 用户链路下行


        }
    }

    return true;
}

void update_link(int time, const vector<string>& selected_link, double need_capacity, Info& info) {
    int length = (int)selected_link.size();
    for (size_t i = 0; i < selected_link.size(); ++i) {
        if (i == 0) {
            // 用户-卫星 上行
            if (facility_name.find(selected_link[length - 1]) != facility_name.end()) {
                info.satellite[time][selected_link[i]].remaining_user_up_capacity -= need_capacity;
            }
            else {
                // 地面站-卫星 上行
                info.satellite[time][selected_link[i]].remaining_faci_up_capacity -= need_capacity;
            }
        }
        else {
            if (info.users.find(selected_link[i]) == info.users.end()) {
                // 星间链路
                if (facility_name.find(selected_link[i]) == facility_name.end() &&
                    facility_name.find(selected_link[i - 1]) == facility_name.end()) {
                    info.satellite_link[time].at(selected_link[i - 1])[selected_link[i]].remain -= need_capacity;
                }
                // 馈电下行
                else if (facility_name.find(selected_link[i - 1]) == facility_name.end() &&
                    facility_name.find(selected_link[i]) != facility_name.end()) {
                    double distance = Distance3D(facility_position.at(selected_link[i]).p,info.satellite_position[time].at(selected_link[i - 1]).p);
                    double need_power = calculate_sat_faci_need_power(need_capacity,sat_antenna.FeedLink_DL_Bandwith,distance);
                    info.satellite[time].at(selected_link[i - 1]).remaining_faci_down_capacity -= need_capacity;
                    info.satellite[time].at(selected_link[i - 1]).remaining_faci_power -= need_power;
                }
                // 馈电上行
                else if (facility_name.find(selected_link[i - 1]) != facility_name.end() &&
                    facility_name.find(selected_link[i]) == facility_name.end()) {
                    info.satellite[time][selected_link[i]].remaining_faci_up_capacity -= need_capacity;
                }
            }
            else {
                // 用户链路下行
                double distance = Distance3D(info.user_position[time][selected_link[i]].p,
                    info.satellite_position[time].at(selected_link[i - 1]).p);
                double need_power = calculate_sat_user_need_power(need_capacity, sat_antenna.UsrLink_DL_Bandwith,distance);
                info.satellite[time].at(selected_link[i - 1]).remaining_user_down_capacity -= need_capacity;
                info.satellite[time].at(selected_link[i - 1]).remaining_user_power -= need_power;
            }
        }
    }
}



const int packet_length = 2048;

// Q-function (similar to Python's erf)
double Q_function(double x) {
    return 0.5 * (1 - erf(x / sqrt(2)));
}

// Calculate packet loss rate from BER and packet length
double calculate_packet_loss(double ber, int packet_length) {
    return 1 - pow(1 - ber, packet_length);
}

// Calculate transmission rate based on SNR and bandwidth
double calculate_rate(double snr, double bandwidth) {
    return CNR2FrequenceEff(snr) * bandwidth;  // Assuming you have a CNR2FrequenceEff function
}

// Calculate required bandwidth based on target rate and SNR
double calculate_required_bandwidth(double rate, double snr) {
    return rate / log2(1 + snr);
}

// Calculate packet loss based on SNR
double calculate_loss(double SNR) {
    double BER = Q_function(sqrt(2 * SNR));
    double loss = 1 - pow(1 - BER, packet_length);

    // Ensure loss is not negative
    return max(0.0, round(loss * 1000000) / 1000000);  // rounding to 6 decimal places
}
// 计算用户上行链路SNR
double calculate_user_up_SNR(double distance) {
    // FSPL
    double Lp = 20 * log10(4 * global_PI * distance * sat_antenna.UsrLink_UL_Frequency / global_C_light);
    double snr = user_antenna.EIRP + sat_antenna.UsrLink_UL_G_t - global_K - 10 * log10(sat_antenna.UsrLink_UL_Bandwith) - Lp;
    return snr;
}

// 计算馈电链路上行SNR
double calculate_faci_up_SNR(double distance) {
    double Lp = 20 * log10(4 * global_PI * distance * sat_antenna.FeedLink_UL_Frequency / global_C_light);
    // fac_antenna.Pt 单位是W
    double Pt_dBW = 10 * log10(fac_antenna.Pt);

    double eirp = Pt_dBW + fac_antenna.G_t - Lp;
    double snr = eirp + sat_antenna.FeedLink_UL_G_t - global_K - 10 * log10(sat_antenna.FeedLink_UL_Bandwith);
    return snr;
}

// 计算卫星用户链路下行SNR
double calculate_user_down_SNR(double distance, double power) {
    // FSPL
    double Lp = 20 * log10(4 * global_PI * distance * sat_antenna.UsrLink_DL_Frequency / global_C_light);
    // power 单位是W
    double Pt_dBW = 10 * log10(power);
    // (dBW)EIRP = Pt - Lp + Gt
    double eirp = Pt_dBW + sat_antenna.UserLink_DL_Gt - Lp;
    double g_t = user_antenna.G_r - 10 * log10(user_antenna.T) - user_antenna.NF;
    // (dBW)snr = EIRP + G/T - K - B
    double snr = eirp + g_t - global_K - 10 * log10(sat_antenna.UsrLink_DL_Bandwith);
    return snr;
}

// 计算卫星馈电链路下行SNR
double calculate_faci_down_SNR(double distance, double power) {
    double Lp = 20 * log10(4 * global_PI * distance * sat_antenna.FeedLink_UL_Frequency / global_C_light);
    // power 单位是W
    double Pt_dBW = 10 * log10(power);

    double eirp = Pt_dBW + sat_antenna.FeedLink_DL_Gt - Lp;
    double g_t = fac_antenna.G_r - 10 * log10(fac_antenna.T) - fac_antenna.NF;
    double snr = eirp + g_t - global_K - 10 * log10(sat_antenna.FeedLink_DL_Bandwith);
    return snr;
}

// 计算馈电链路所需的功率
double calculate_sat_faci_need_power(double rate, double bandwidth, double distance) {
    double Lp = 20 * log10(4 * global_PI * distance * sat_antenna.FeedLink_DL_Frequency / global_C_light);
    double snr = pow(2, rate / bandwidth) - 1;
    double g_t = fac_antenna.G_r - 10 * log10(fac_antenna.T) - fac_antenna.NF;
    double eirp = snr - g_t + global_K + 10 * log10(sat_antenna.FeedLink_DL_Bandwith);
    double Pt = eirp + Lp - sat_antenna.FeedLink_DL_Gt;
    double Pt_W = pow(10, Pt / 10);
    return Pt_W;
}

// 计算用户链路所需的功率
double calculate_sat_user_need_power(double rate, double bandwidth, double distance) {
    //       G/T = Gr - NF - 10log10(T)
    //  (dBW)snr = EIRP + G/T - K - B
    //  (dBW)EIRP = Pt - Lp + Gt
    double Lp = 20 * log10(4 * global_PI * distance * sat_antenna.UsrLink_DL_Frequency / global_C_light);
    double snr = pow(2, rate / bandwidth) - 1;
    double g_t = user_antenna.G_r - 10 * log10(user_antenna.T) - user_antenna.NF;
    double eirp = snr - g_t + global_K + 10 * log10(sat_antenna.UsrLink_DL_Bandwith);
    double Pt = eirp + Lp - sat_antenna.UserLink_DL_Gt;
    double Pt_W = pow(10, Pt / 10);
    return Pt_W;
}

string joinWithCommas(const vector<string>& strings) {
    if (strings.empty()) {
        return "";
    }

    ostringstream oss;
    for (size_t i = 0; i < strings.size(); ++i) {
        oss << strings[i];
        if (i != strings.size() - 1) {
            oss << ",";
        }
    }
    return oss.str();
}