//
// Created by TXW on 24-12-15.
//
#include "dijkstra.h"




void Dijskra::initGraph(const vector<vector<double>>& weight)
{
	if (init_flag)
	{
		VexList.clear();
		pathRecord.clear();
		distRecord.clear();
		doneRecord.clear();
	}
	EdgeNode* e = new EdgeNode();
	Vertexs = weight.size();
	for (int i = 0; i < Vertexs; ++i) { //初始化顶点
		VertexList tmp;
		tmp.firstedge = NULL;
		VexList.emplace_back(tmp);
	}
	for (int src = 0; src < Vertexs; ++src) //初始化边
	{
		for (int dst = 0; dst < Vertexs; ++dst)
		{
			if (weight[src][dst] >= INFINITY || src == dst) continue;
			if (VexList[src].firstedge == NULL) {//当前顶点i后面没有顶点
				e = new EdgeNode;
				e->adjvex = dst;
				e->val = weight[src][dst];
				e->next = NULL;
				VexList[src].firstedge = e;
			}
			else {  //当前i后面有顶点
				EdgeNode* p = VexList[src].firstedge;
				while (p->next) {
					p = p->next;
				}
				e = new EdgeNode;
				e->adjvex = dst;
				e->val = weight[src][dst];
				e->next = NULL;
				p->next = e;
			}
		}
	}
	pathRecord.resize(Vertexs);
	distRecord.resize(Vertexs);
	init_flag = true;
}


double Dijskra::dijskra(int src, int dst)
{
	if (!init_flag)return INFINITY;
	if (doneRecord.count(src) != 0)return distRecord[src][dst];
	vector<double> dist(Vertexs, INFINITY);
	vector<int> path(Vertexs);
	for (int i = 0; i < Vertexs; ++i)path[i] = i;
	// 基于优先队列实现堆，根据边的权值进行排序
	priority_queue<EdgeNode> min_edge;
	vector<int> vis(Vertexs, 0);
	EdgeNode* tmp = VexList[src].firstedge;
	while (tmp)
	{
		dist[tmp->adjvex] = tmp->val;
		path[tmp->adjvex] = src;
		min_edge.push(EdgeNode(tmp->adjvex, tmp->val));
		tmp = tmp->next;
	}
	dist[src] = 0;
	vis[src] = 1;
	while (!min_edge.empty())
	{
		EdgeNode edge = min_edge.top();
		min_edge.pop();
		int cur_node = edge.adjvex;
		if (vis[cur_node] == 1)continue;
		vis[cur_node] = 1;
		EdgeNode* tmp_edge = VexList[cur_node].firstedge;
		while (tmp_edge)
		{
			if (dist[tmp_edge->adjvex] > dist[cur_node] + tmp_edge->val)
			{
				dist[tmp_edge->adjvex] = dist[cur_node] + tmp_edge->val;
				path[tmp_edge->adjvex] = cur_node;
				if (vis[tmp_edge->adjvex] == 0)
				{
					min_edge.push(EdgeNode(tmp_edge->adjvex, dist[tmp_edge->adjvex]));
				}
			}
			tmp_edge = tmp_edge->next;
		}
	}
	distRecord[src] = dist;
	pathRecord[src] = path;
	doneRecord.insert(make_pair(src, 0));
	return dist[dst];
}
void Dijskra::DFS(vector<int>& path, int src, int dst)
{
	if (src == dst)
	{
		path.emplace_back(src);
		return;
	}
	DFS(path, src, pathRecord[src][dst]);
	path.emplace_back(dst);
}
void Dijskra::findWaypoint(vector<int>& path, int src_node, int dst_node)
{
	if (!init_flag || doneRecord.count(src_node) == 0)return;
	DFS(path, src_node, dst_node);
}