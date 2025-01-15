#include "max.h"
#include "qoe.h"
#include "resource.h"
#include <filesystem>

using namespace std;

int main()
{
    int choice;
    cin >> choice;
    filesystem::path current_path = filesystem::current_path();
    cout << "当前工作目录: " << current_path << endl;
    if (choice == 1)
    {
        int len,choice,plus;
        cin>>len>>choice>>plus;
        calculate_all_qoe(current_path,len,choice,plus);
    }
    else if (choice == 2)
    {
        int len,choice,plus;
        cin>>len>>choice>>plus;
        calculate_all_max(current_path,len,choice,plus);
    }else if(choice == 3){
        int num;
        cin>>num;
        Info info;
        info.init(2000,-1,current_path);
        info.calculate_all_link(num);
    }
    return 0;
}