// worm.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
using namespace cv; 
using namespace std;

const int COLOR_RANGE = 256;
const int BLANK = 255;
const int BLACK = 0;
const double INF = 1e100;
const double MIN_WORM_SIZE_RATE = 0.1;
const double MIN_WORM_SEG_RATE = 0.2;

typedef vector<int> vi;
typedef vector<vi> vii;

typedef vector<uchar> vu;
typedef vector<vu> vuu;

inline double sqr(double x)
{
	return x * x;
}

const int totalDir = 8;
const int dir[totalDir][2] = {{-1, -1},{-1, 0}, {-1,1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}};

template <class T>
inline void initVec(vector <vector <T> > &x, int n, int m, T v)
{
	x.resize(n);
	for (int i = 0; i < n; i++) x[i].resize(m, v);
}

class Worm
{
public:
	Worm(int s): size(s)
	{
	}
	vector <Point> margin;
	int size;
};

class Edge
{
public:
	Edge(int x, int y, int length, bool alive, int father, bool xused, bool yused) 
		:x(x), y(y), length(length), alive(alive), father(father), xused(xused), yused(yused)
	{
	}
	Edge(int x, int y, int length) :x(x), y(y), length(length), alive(true), father(-1), xused(false), yused(false)
	{
	}
	int length, x, y;
	bool xused, yused;
	int father;
	bool alive;
};

typedef set <int> Edgelist;

class Graph
{
public:
	int n,m;
	vector <Edgelist> eList;
	vector <Edge> edges;
	vector <Point> nodes;
	bool findEdge(int x, int y, Edge* e) {
		for (vector  <Edge>::iterator it = edges.begin(); it != edges.end(); it++)
			if ( ((it->x == x) && (it->y == y)) || 
				((it->x == y) && (it->y == x)) ) {
					e = &(*it);
					return true;
			}
		return false;
	}
	int addNode(Point &p)
	{
		nodes.push_back(p);
		n++;
		eList.push_back(Edgelist());
		return nodes.size() - 1;
	}
	int addNode(int x, int y)
	{
		nodes.push_back(Point(x, y));
		n++;
		eList.push_back(Edgelist());
		return nodes.size() - 1;
	}
	int addEdge(int x, int y, int length)
	{
		edges.push_back(Edge(x,y,length));
		int num = edges.size() - 1;
		eList[x].insert(num);
		eList[y].insert(num);
		m++;
		return num;
	}
	int addEdge(Edge e)
	{
		edges.push_back(e);
		int num = edges.size() - 1;
		eList[e.x].insert(num);
		eList[e.y].insert(num);
		m++;
		return num;
	}
	void falselyRemoveEdge(Edge &tedge)
	{
		if (tedge.alive)
		{
			tedge.alive = false;
			m--;
		}
	}
	void unionNode(int x, int y)
	{
		printf("%d %d\n", x, y);
		eList[x].insert(eList[y].begin(), eList[y].end());
		for (vector<Edge>::iterator it = edges.begin(); it != edges.end(); it++)
		{
			if (it -> x == y)
			{
				it -> x = x;
				if (it -> y == x) falselyRemoveEdge(*it);
			}
			if (it -> y == y)
			{
				it -> y = x;
				if (it -> x == x) falselyRemoveEdge(*it);
			}
		}
		n--;
	}
	Graph(): n(0), m(0)
	{
		edges.clear();
		eList.clear();
		nodes.clear();
	}
	void clear()
	{
		n = 0;
		m = 0;
		eList.clear();
		edges.clear();
		nodes.clear();
	}
	void cleanEdge(int lowest)
	{
		for (vector <Edge>::iterator it = edges.begin(); it != edges.end(); it++)
		{
			if ( it -> alive && it-> length < lowest)
			{

				unionNode(it->x, it -> y);
				falselyRemoveEdge(*it);
			}
		}
	}
	void printInfo()
	{
		printf("n = %d m = %d", n, m);
		printf("nodes:\n");
		for (vector <Point>::iterator it = nodes.begin(); it != nodes.end(); it++)
			printf("x:%d ,y:%d\n", it -> x, it -> y);
		printf("edges:\n");
		for (vector <Edge>::iterator it = edges.begin(); it != edges.end(); it++)
			if (it -> alive)
				printf("x:%d ,y:%d\n", it -> x, it -> y);
		printf("List:\n");
		for (int i = 0; i < eList.size(); i++)
		{
			printf("node %d:", i);
			for (Edgelist::iterator it = eList[i].begin(); it != eList[i].end(); it++)
				if (edges[*it].alive)
					printf("%d ", *it);
			printf("\n");
		}
	}
private:

};

class WormMat :public Mat
{

public:
	vector <int> region;
	vector <Worm> worms;
	vector< vector <uchar>> sign;
	int countNumber()
	{
		Graph g;
		vii visit;
		initVec(visit, rows, cols, 0);
		int ret = 0;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (!visit[i][j] && this -> at<uchar>(i,j) == BLACK)
				{
					g.clear();
					int maxLength = 0;

					constructGraph(*this, g, visit, i, j, Edge(-1,-1,0), maxLength, false);
					g.cleanEdge(maxLength * MIN_WORM_SEG_RATE);
					g.printInfo();
					int tmp_number = countNumber(g);
					printf("Worm number: %d---------\n", tmp_number);
					ret += tmp_number;
				}
		return ret;
	}

	WormMat () : Mat()
	{
		regionIsSet = 0;
		this->sign.clear();
	}
	WormMat (const Mat &b) : Mat(b)
	{
		regionIsSet = 0;
		sign.clear();
	}
	~WormMat()
	{
		~Mat();
	}
	WormMat extract(int layer, int s, int e)
	{
		s--;
		e--;
		if (regionIsSet != layer) countRegion(layer);
		WormMat ret = *this;
		for (Mat_<uchar>::iterator it = ret.begin<uchar>(); it!= ret.end<uchar>(); ++it)
		{
			int t = 0;
			while (*it > region[t + 1]) t++;
			if (t < s || t > e) *it = 255;
		}
		return ret;
	}

	void countRegion(int layer)
	{
		regionIsSet = layer;
		int color[COLOR_RANGE], tot[COLOR_RANGE];
		vii pre;
		vector <vector <double> >opt;
		region.resize(layer + 1, 0);
		initVec(pre, COLOR_RANGE, layer + 1, -1);
		initVec(opt, COLOR_RANGE, layer + 1, 0.0);
		double sumx[COLOR_RANGE], sumx2[COLOR_RANGE];
		memset(color, 0, sizeof color);

		for (Mat_<uchar>::iterator it= this -> begin<uchar>(); it!= this -> end<uchar>(); ++it)
			color[*it]++;
		/*for (int i = 0; i < COLOR_RANGE; i++)
		{
			printf("%d %d\n", i, color[i]);
		}*/
		pretreatment(tot, sumx, sumx2, color);
		for (int i = 1; i < COLOR_RANGE; i++)
			opt[i][1] = tot[i] ? (sumx2[i] / tot[i] - sqr(sumx[i] / tot[i])) : 0;
		for (int i = 1; i < COLOR_RANGE; i++)
			for (int k = 2; k <= layer; k++) opt[i][k] = INF;
		for (int i = 2; i < COLOR_RANGE; i++)
			for (int j = 1; j < i; j++)
			{
				double variance = 0;
				int n = tot[i] - tot[j - 1];
				if (n > 0)
				{
					double ex2 = (sumx2[i] - sumx2[j - 1]) / n, e2x = sqr((sumx[i] - sumx[j - 1]) / n);
					variance = ex2 - e2x;
				}
				for (int k = 2; k <= layer; k++)
				{
					if (opt[j][k - 1] + variance < opt[i][k])
					{
						opt[i][k] = opt[j][k - 1] + variance;
						pre[i][k] = j;
					}
				}
			}
			//printf("%d\n",pre[COLOR_RANGE - 1][layer]);
			for (int p = COLOR_RANGE - 1, i = layer; i > 1; i--, p = pre[p][i])
				region[i - 1] = pre[p][i];
			region[0] = 0;
			region[layer] = COLOR_RANGE - 1;
			for (int i = 0; i <= layer; i++) printf("%d ",region[i]);
	}

	void thinning()
	{
		static const int erasetable[256]={
0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0
		};
		bool Finished = false;
		while(!Finished)
		{
			Finished = true;
			for (int y = 1; y < cols - 1; y++)
			{ 
				int x = 1; 
				while(x< rows - 1)
				{
					if(this -> at<uchar>(x,y) == BLACK)
					{
						if( (this -> at<uchar>(x - 1,y)==BLANK)|| (this -> at<uchar>(x + 1,y)==BLANK))
						{
							int num = 0;
							for (int i = 0; i < totalDir; i++)
								num |= ((this -> at<uchar>(x + dir[i][0], y + dir[i][1])) / 255) << i;
							if(erasetable[num] == 1)
							{
								this -> at<uchar>(x,y) = BLANK;
								Finished = false;
								x++;
							}
						}
					}
					x++;
				}
			}
			for (int x = 1;x< rows - 1; x++)
			{ 
				int y = 1;
				while(y < cols - 1)
				{
					if(this -> at<uchar>(x,y)== BLACK)
					{
						if( (this -> at<uchar>(x,y - 1) == BLANK)|| (this -> at<uchar>(x,y + 1) == BLANK))
						{
							int num = 0;
							for (int i = 0; i < totalDir; i++)
								num |= ((this -> at<uchar>(x + dir[i][0], y + dir[i][1])) / 255) << i;
							if(erasetable[num] == 1)
							{
								this -> at<uchar>(x,y) = BLANK;
								Finished = false;
								y++;
							}
						}
					}
					y++;
				}
			} 
		}
	}

	void split(WormMat &afterClean, WormMat &margin)
	{
		Mat &a = *this;
		int n = a.rows, m = a.cols;
		initVec(sign ,n,m, uchar(0));
		int cnt = 0, mSize = 0, tar = 0; 
		vector <int> size;
		size.clear();
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				if (!sign[i][j] && a.at<uchar>(i,j) < BLANK)
				{
					int nowSize = 0;
					myFloodFill(a, i, j, ++cnt, sign, nowSize);
					size.push_back(nowSize);
					if (nowSize > mSize)
					{
						mSize = nowSize;
						tar = cnt;
					}
				}
		cnt = 0;
		vii isMargin;
		initVec(isMargin, rows, cols, 0);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				if (sign[i][j] > cnt)
				{
					//printf("%d %d %d\n",size[sign[i][j] - 1], mSize, sign[i][j]);
					if (size[sign[i][j] - 1] < mSize * MIN_WORM_SIZE_RATE) fillInt(sign, i, j, 0);
					else
					{
						size[++cnt - 1] = size[sign[i][j] - 1];
						fillInt(sign, i, j, cnt);
						findMargin(sign, i, j, isMargin, 0);
					}
				}
		afterClean = Mat(rows, cols, CV_8UC1);
		margin = Mat(rows, cols, CV_8UC1);
		Mat_<uchar>::iterator it= afterClean.begin<uchar>();
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (sign[i][j]) *it++ = BLACK;
				else *it++ = BLANK;
		it= margin.begin<uchar>();
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (isMargin[i][j])
					*it++ = BLACK; 
				else *it++ = BLANK;
	}

private:
	int regionIsSet;

	void constructGraph(Mat &a, Graph &g, vii &visit, int x, int y, Edge e, int &maxLength, bool newNode)
	{
		e.length++;
		maxLength = max(maxLength, e.length);
		int cnt = 0;
		for (int i = 0; i < totalDir; i++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx <rows && ty >= 0 && ty < cols && a.at<uchar>(tx,ty) == BLACK)
			{
				if (!newNode && visit[tx][ty] > 1)
				{
					e.y = visit[tx][ty] - 2;
					g.addEdge(e);
					return;
					//return;
				}
				else
					if (!visit[tx][ty]) cnt++;
			}
		}
		newNode = cnt == 0 || cnt > 1 || e.x == -1;
		if (newNode)
		{
			e.y = g.addNode(x,y);
			//e.length = dist(g)
			if (e.x != -1) g.addEdge(e);
			visit[x][y] = e.y + 2;
			if (cnt == 0) return;
			e.x = e.y;
			e.y = -1;
			e.length = 0;
		}
		else
			visit[x][y] = 1;
		for (int i = 0; i < totalDir; i++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx <rows && ty >= 0 && ty < cols && a.at<uchar>(tx,ty) == BLACK && !visit[tx][ty]) constructGraph(a, g, visit, tx, ty, e, maxLength, newNode);
		}
	}

	int countNumber(Graph &g)
	{
		while (true) {
			double opt_cost = 0; // must be an acute angle if the combination is triggered
			// find the next edge to expand
			Edge* first_edge = NULL, * second_edge = NULL;
			for (vector <Edge>::iterator it = g.edges.begin(); it != g.edges.end(); it++) {
				if (!it->xused) 
					testProlongation(g, &(*it), it->x, opt_cost, first_edge, second_edge);
				if (!it->yused) 
					testProlongation(g, &(*it), it->y, opt_cost, first_edge, second_edge);
			}
			printf("cost:%lf first:%d second:%d\n", opt_cost, first_edge == NULL, second_edge == NULL);
			if (first_edge != NULL)
				printf("first: %d %d %d %d\n", first_edge->x, first_edge->y, first_edge->xused, first_edge->yused);
			if (second_edge!= NULL)
				printf("second: %d %d %d %d\n", second_edge->x, second_edge->y, second_edge->xused, second_edge->yused);

			// if no edge can be extend, return the counted value
			if (first_edge == NULL) {
				int count = 0;
				for (vector <Edge>::iterator it = g.edges.begin(); it != g.edges.end(); it++)
					count += it->alive;
				return count;
			}
			/*
				int count = 0;
				for (vector <Edge>::iterator it = g.edges.begin(); it != g.edges.end(); it++)
					count += it->alive;
				printf("count: %d\n", count);
				*/
			// else combine the edges
			int common_node = getCommonNode(*first_edge, *second_edge);
			int first_other_node = first_edge->x + first_edge->y - common_node;
			int second_other_node = second_edge->x + second_edge->y - common_node;
			// decide which edge to delete (the one with a node with one degree, else select the one with ...?)
			//TODO here, for the second condition, we simply select the second node
			Edge* edge_to_del = NULL;
			Edge* edge_remained = NULL;
			//printf("common node %d\n", common_node);
			printf("edgesize of %d:%d used info %d\n",first_other_node, getEdgelistSize(g, first_other_node), getUsedInfo(*first_edge, common_node));
			if (getEdgelistSize(g, first_other_node) == 1 && !getUsedInfo(*first_edge, common_node)) {
				edge_to_del = first_edge;
				edge_remained = second_edge;
			}
			else {
				edge_to_del = second_edge;
				edge_remained = first_edge;
			}
			if (edge_to_del != NULL)
				printf("to_del: %d %d %d %d\n", edge_to_del->x, edge_to_del->y, edge_to_del->xused, edge_to_del->yused);
			if (edge_remained!= NULL)
				printf("remained: %d %d %d %d\n", edge_remained->x, edge_remained->y, edge_remained->xused, edge_remained->yused);

			// this may not be necessary to assign ``used'' information on the edge to delete
			if (edge_to_del->x == common_node)
				edge_to_del->xused = true;
			else
				edge_to_del->yused = true;
			edge_to_del->alive = false;

			if (edge_remained->x == common_node){
				if (edge_remained->xused == true) {// x has been combined
					//printf("Get here\n");
					// copy edge_remained, and create a new node
					int new_node_index = g.addNode(g.nodes[common_node]);
					Edge new_edge(*edge_remained);
					new_edge.x = new_node_index;
					//TODO fix up the chaining information of edges
					// new_edge.father = (index of) edge_to_del
					g.addEdge(new_edge);
				}
				edge_remained->xused = true; // x has not been combined previously
			} else if (edge_remained->y == common_node) {
				if (edge_remained->yused == true) { // y has been combined
					//printf("Get here\n");
					int new_node_index = g.addNode(g.nodes[common_node]);
					Edge new_edge(*edge_remained);
					new_edge.y = new_node_index;
					//TODO fix up the chaining information of edges
					// new_edge.father = (index of) edge_to_del
					g.addEdge(new_edge);
				}
				edge_remained->yused = true;
			}
		}
		return 0;
	}

	bool getUsedInfo(const Edge& e, const int node_index) {
		//printf("In used info, node_index:%d, e.x:%d, x_used%d; e.y:%d, y_used:%d\n", node_index, e.x, e.xused, e.y, e.yused);
		if (e.x == node_index)
			return e.xused;
		else
			return e.yused;
	}

	void putUsedInfo(Edge& e, const int node_index) {
		if (e.x == node_index)
			e.xused = true;
		else if (e.y == node_index)
			e.yused = true;
		else
			printf("This edge does not contain the specified node.\n");
			exit(0);
	}

	// try prolonging the edge from ``edge_node''
	bool testProlongation(Graph& g, Edge* it, int edge_node, double &opt_cost, Edge *&first_edge, Edge *&second_edge) {
		bool modified = false;
		for (Edgelist::iterator other_edge_index = g.eList[edge_node].begin(); 
			other_edge_index != g.eList[edge_node].end();
			other_edge_index++) {
				// if the edge does not exist or if the other edge is itself, then continue
				if (!g.edges[*other_edge_index].alive ||
					&g.edges[*other_edge_index] == it)
					continue;
				// if the two edges have combined contacted nodes
				if (getUsedInfo(*it, edge_node) && getUsedInfo(g.edges[*other_edge_index], edge_node))
					continue;
				double tmp_cost = calcCost(g, *it, g.edges[*other_edge_index]);
				if (tmp_cost < opt_cost) {
					first_edge = &(*it);
					second_edge = &g.edges[*other_edge_index];
					opt_cost = tmp_cost;
					modified = true;
				}
		}
		return modified;
	}

	int getCommonNode(const Edge& first_edge, const Edge& second_edge) {
		int common_node = -1;
		if ((first_edge.x == second_edge.x) || (first_edge.x == second_edge.y))
			common_node = first_edge.x;
		else
			common_node = first_edge.y;
		return common_node;
	}

	double calcCost(Graph& g, const Edge& first_edge, const Edge& second_edge) {
		// should consider the previous combining history and edge's local information (derivative at a small distance)
		int common_node = getCommonNode(first_edge, second_edge);
		int first_other_node = first_edge.x + first_edge.y - common_node;
		int second_other_node = second_edge.x + second_edge.y - common_node;
		// calculate the cosine value of angle inbetween
		double x1 = g.nodes[common_node].x - g.nodes[first_other_node].x ,
			y1 = g.nodes[common_node].y - g.nodes[first_other_node].y;
		double x2 = g.nodes[common_node].x - g.nodes[second_other_node].x ,
			y2 = g.nodes[common_node].y - g.nodes[second_other_node].y;
		double retVal = (x1*x2 + y1*y2) / sqrt(x1*x1+y1*y1) / sqrt(x2*x2+y2*y2);
		// if one edge is of degree zero, then it should be considered first (let a zero calcCost if the angle is of good value)
		// TODO may need to modify???
		if (!getEdgelistSize(g, first_other_node) || !getEdgelistSize(g, second_other_node))
			return -2; // optimal
		else
			return retVal;
	}

	int getEdgelistSize(Graph& g, int index) {
		int count = 0;
		for (Edgelist::iterator it = g.eList[index].begin(); it != g.eList[index].end(); it++)
			count += g.edges[*it].alive;
		return count;
	}

	void pretreatment(int *tot, double *sumx, double *sumx2, int *color)
	{
		memset(tot, 0, COLOR_RANGE * sizeof(int));
		memset(sumx, 0, COLOR_RANGE * sizeof(int));
		memset(sumx2, 0, COLOR_RANGE * sizeof(int));
		tot[0] = color[0];
		sumx[0] = 0;
		sumx2[0] = 0;
		for (int i = 1; i < COLOR_RANGE; i++)
		{
			tot[i] = tot[i - 1] + color[i];
			sumx[i] = sumx[i - 1] + double(color[i]) * i;
			sumx2[i] = sumx2[i - 1] + double(color[i]) * i * i;
		}
	}

	void myFloodFill(Mat &a, int x, int y, int cnt, vuu &sign, int &size)
	{
		queue <Point> q;
		q.push(Point(x,y));
		sign[x][y] = cnt;
		while (!q.empty())
		{
			//printf("%d %d %d\n",x,y,cnt);
			x = q.front().x;
			y = q.front().y;
			q.pop();
			size++;
			for (int i = 0; i < totalDir; i++)
			{
				int tx = x + dir[i][0];
				int ty = y + dir[i][1];
				if (tx >= 0 && tx <a.rows && ty >= 0 && ty < a.cols && sign[tx][ty] == 0 && a.at<uchar>(tx,ty) < BLANK)
				{
					sign[tx][ty] = cnt;
					q.push(Point(tx, ty));
				}
			}
		}
		//printf("%d %d\n", a.at<uchar>(x,y));
	}

	void fillInt(vuu &sign, int x, int y, int cnt)
	{
		queue <Point> q;
		q.push(Point(x,y));
		sign[x][y] = cnt;
		while (!q.empty())
		{
			//printf("%d %d %d\n",x,y,cnt);
			x = q.front().x;
			y = q.front().y;
			q.pop();
			for (int i = 0; i < totalDir; i++)
			{
				int tx = x + dir[i][0];
				int ty = y + dir[i][1];
				if (tx >= 0 && tx <rows && ty >= 0 && ty < cols && sign[tx][ty] != cnt && sign[tx][ty])
				{
					sign[tx][ty] = cnt;
					q.push(Point(tx, ty));
				}
			}
		}
	}

	void findMargin(vuu &sign, int x, int y, vii &isMargin, int nowdir)
	{
		//printf("%d %d %d\n", x, y, isMargin[x][y]);
		for (int i = nowdir, j = 0; j < totalDir; j++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx < rows && ty >= 0 && ty < cols && sign[tx][ty] == sign[x][y])
			{
				if (!(isMargin[x][y] & (1 << i)))
				{
					isMargin[x][y] |= 1 << i;
					findMargin(sign, tx, ty, isMargin, (i + totalDir - 1) % totalDir);
				}
				break;
			}
			i = (i + 1) % totalDir;
		}
	}
};

Mat orgImg, imgGray;

void init(String fn)
{
	orgImg = imread(fn); 
	if (orgImg.empty()) 
	{
		fprintf(stderr,"Error: load image failed."); 
		exit(-1); 
	} 
	cvtColor(orgImg, imgGray, CV_RGB2GRAY);
	printf("%d\n", imgGray.type());

}

int _tmain(int argc, _TCHAR* argv[])
{
	init("1_2.png");
	WormMat	img = imgGray;

	imwrite("gray.bmp",img);
	img.countRegion(6);
	medianBlur(img, img, 27);
	imwrite("blur.bmp",img);
	WormMat newimg = img.extract(6,1,3);
	namedWindow("image", CV_WINDOW_AUTOSIZE);
	namedWindow("image2", CV_WINDOW_AUTOSIZE);
	namedWindow("image3", CV_WINDOW_AUTOSIZE);
	imshow("image2", newimg);
	imwrite("extract.bmp",newimg);
	WormMat afterClean, margin;
	newimg.split(afterClean, margin);
	imwrite("clean.bmp",afterClean);

	imshow("image", afterClean);
	afterClean.thinning();

	imwrite("thinning.bmp",afterClean);
	afterClean.countNumber();
	imshow("image3", afterClean);
	//imwrite("test2.bmp",afterClean);
	imshow("image2", margin);
	waitKey(); 
	return 0;
}

