//
// Created by TXW on 24-12-15.
//

#ifndef DIJKSTRA_H
#define DIJKSTRA_H
#pragma once
#include <iostream>
#include <cstdio>
#include <cstring>
#include <queue>
#include <vector>
#include <unordered_map>
#include <cmath>


using namespace std;

// 定义边节点
typedef struct EdgeNode {
    int adjvex; // 顶点下标
    struct  EdgeNode* next; // 下一条边的指针
    double val;    // 当前边的代价

    EdgeNode() {
        adjvex = -1;
        val = 0;
        next = nullptr;
    };
    // 构造函数
    EdgeNode(int _adjvex, double _val)
    {
        adjvex = _adjvex;
        val = _val;
        next = nullptr;
    }
    // 重载<操作符
    bool operator<(const EdgeNode& other) const
    {
        // 对小于运算符进行重载，如果两个值相等，那么继续判断point，如果不等，则返回
        if (this->val == other.val)
        {
            return this->adjvex < other.adjvex;
        }
        else
        {
            return this->val > other.val;
        }
    }

} EdgeNode;


// 定义顶点表
struct VertexList
{
    EdgeNode* firstedge;  //用来存储当前顶点的下一个顶点

    VertexList() {
        firstedge = nullptr;
    };
    ~VertexList() {
        EdgeNode* tmp = firstedge;
        while (firstedge)
        {
            firstedge = firstedge->next;
            delete tmp;
            tmp = firstedge;
        }
    };
};


class Dijskra
{
public:
    Dijskra() :Vertexs(0), init_flag(false) {};
    void initGraph(const vector<vector<double>>& weight);
    double dijskra(int src_node, int dst);
    const vector<double>& getDist(int src) { return distRecord[src]; };
    const vector<int>& getPrePath(int src) { return pathRecord[src]; };
    void findWaypoint(vector<int>& path, int src_node, int dst_node);
private:
    bool init_flag;	//记录是否初始化过
    int Vertexs;	//顶点数量
    vector<VertexList> VexList;				//顶点集合
    vector<vector<int>>	pathRecord;			//记录已经计算过的最短路径
    vector<vector<double>> distRecord;		//记录已经计算过的最短距离
    unordered_map<int, int> doneRecord;		//记录已寻路过的节点

    void DFS(vector<int>& path, int src, int dst);
};


#endif //DIJKSTRA_H
