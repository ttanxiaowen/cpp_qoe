#pragma once


#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "resource.h"


double play_time_comparison(double x0, double x1);
double exp_loss(double x);
double exp_time(double x);
double custom_function(double x);
double calculate_qoe(double rate, double need_rate, double delay, double loss);


double calculate_user_up_SNR(double distance);
double calculate_loss(double snr);
double calculate_rate(double snr, double bandwidth);
double calculate_sat_faci_need_power(double rate, double bandwidth, double distance);
double calculate_user_down_SNR(double distance, double power);
double calculate_sat_user_need_power(double rate, double bandwidth, double distance);




void update_link(int time, const vector<string>& selected_link, double need_capacity, Info& info);
bool can_meet_link_requirements(int time, int index, const vector<string>& selected_link, Info& info, double* rate, double* loss, vector<int>& unsatis,int episode);
string joinWithCommas(const vector<string>& strings);