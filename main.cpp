// //
// // Created by TXW on 24-12-15.
// //
// #include <iostream>
// #include <filesystem>
//
// std::filesystem::path find_root_directory(const std::filesystem::path& start_path) {
//     std::filesystem::path current_path = start_path;
//
//     while (current_path != current_path.root_path()) {
//         // 假设项目根目录包含名为 CMakeLists.txt 的文件
//         if (std::filesystem::exists(current_path / "CMakeLists.txt")) {
//             return current_path;
//         }
//         current_path = current_path.parent_path();  // 返回上级目录
//     }
//
//     return current_path;  // 如果未找到根目录，返回当前目录
// }
//
// int main() {
//     // 获取当前工作目录
//     std::filesystem::path current_path = std::filesystem::current_path();
//     std::filesystem::path project_root = find_root_directory(current_path);
//
//     std::cout << "123 " << project_root << std::endl;
//
//     return 0;
// }