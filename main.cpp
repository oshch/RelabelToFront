#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>

using namespace std;

const int INF = 1e9;
const long long INFLL = (long long)1e18;

class Vertex {
private:
	int dist;
	int edges_head;
	long long excess;
	int label;
public:
	Vertex() : dist(INF), edges_head(0) {}

	int getDist();
	int getHead();
	void setDist(int x);
	void increaseHead();
	void initHead();

	long long getExcess();
	int getLabel();
	void setExcess(long long x);
	void setLabel(int x);
	void increaseExcess(long long x);
};

class Edge {
protected:
	int from, to;
	int id;
public:
	Edge() {}
	Edge(int _from, int _to, int _id) : from(_from), to(_to), id(_id) {}
	int getTo();
	int getId();
	int getFrom();
};

class FlowEdge : public Edge {
private:
	int flow, capacity;
	FlowEdge* rev;
public:
	FlowEdge() : rev(NULL) {}
	FlowEdge(int _from, int _to, int _id, int _flow, int _capacity, FlowEdge* _rev) : Edge(_from, _to, _id) {
		flow = _flow;
	 	capacity = _capacity;
		rev = _rev;
	}

	void setRev(FlowEdge* _rev);
	void push(int x);
	int getG();
	int getFlow();
};

class FlowGraph {
private:
	int EdgeCount;
	vector<vector<FlowEdge*> > E;
	vector<Vertex*> V;

public:
	FlowGraph() {}

	int getSize();
	void get();

	void printFlow();

	void push(FlowEdge* e);
	void relabel(int a);
	void initPreflow(int S);

	void discharge(int u, int S, int T, queue<int>& q, vector<int>& used);
	long long RelabelToFront(int S, int T);
};

class Graph {
private:
	int EdgeCount;
	vector<vector<Edge*> > E;
	vector<Vertex*> V;

public:
	Graph() {}
	Graph(const vector<vector<FlowEdge*> >& FlowE, const vector<Vertex*>& FlowV, int minG);
	int getSize();
};

int Vertex::getDist() {
	return dist;
}

int Vertex::getHead() {
	return edges_head;
}

void Vertex::setDist(int x) {
	dist = x;
}

void Vertex::increaseHead() {
	edges_head++;
}

void Vertex::initHead() {
	edges_head = 0;
}

long long Vertex::getExcess() {
	return excess;
}

int Vertex::getLabel() {
	return label;
}

void Vertex::setExcess(long long x) {
	excess = x;
}

void Vertex::setLabel(int x) {
	label = x;
}

void Vertex::increaseExcess(long long x) {
	excess += x;
}

int Edge::getTo() {
	return to;
}

int Edge::getId() {
	return id;
}

int Edge::getFrom() {
	return from;
}

int FlowEdge::getFlow() {
	return flow;
}

void FlowEdge::setRev(FlowEdge* _rev) {
	rev = _rev;
}

void FlowEdge::push(int x) {
	this->flow += x;
	this->rev->flow -= x;
}

int FlowEdge::getG() {
	return capacity - flow;
}

int Graph::getSize() {
	return E.size();
}

Graph::Graph(const vector<vector<FlowEdge*> >& FlowE, const vector<Vertex*>& FlowV, int minG) {
	E.resize(FlowE.size());
	EdgeCount = 0;
	for (int i = 0; i < (int)FlowE.size(); i++) {
		E[i].resize(0);
		for (int j = 0; j < (int)FlowE[i].size(); j++) {
			if (FlowE[i][j]->getG() >= minG) {
				int y = FlowE[i][j]->getTo();
				E[i].push_back(new Edge(i, y, FlowE[i][j]->getId()));
				EdgeCount++;
			}
		}
	}

	V = FlowV;
}

int FlowGraph::getSize() {
	return E.size();
}

void FlowGraph::get() {
	int n, m;
	cin >> n >> m;
	EdgeCount = m;
	E.resize(n);
	V.resize(n);
	for (int i = 0; i < n; i++) {
		V[i] = new Vertex();
	}

	for (int i = 0; i < m; i++) {
		int a, b, c;
		cin >> a >> b >> c;
		a--; b--;


		E[a].push_back(new FlowEdge(a, b, i + 1, 0, c, NULL));
		E[b].push_back(new FlowEdge(b, a, 0, 0, 0, E[a][ E[a].size() - 1 ]));
		E[a][ E[a].size() - 1 ]->setRev(E[b][ E[b].size() - 1 ]);
	}
}

void FlowGraph::push(FlowEdge* e) {
 	int a = e->getFrom();
 	int b = e->getTo();
 	long long d = min(V[a]->getExcess(), (long long)e->getG());
 	e->push(d);
 	V[a]->increaseExcess(-d);
 	V[b]->increaseExcess(d);
}

void FlowGraph::relabel(int a) {
	int ans = INF;
	for (int i = 0; i < (int)E[a].size(); i++) {
		if (E[a][i]->getG() > 0) {
			ans = min(ans, V[ E[a][i]->getTo() ]->getLabel() + 1);
		}
	}
	V[a]->setLabel(ans);
}

void FlowGraph::initPreflow(int S) {
	for (int i = 0; i < (int)V.size(); i++) {
		V[i]->setExcess(0);
		V[i]->setLabel(0);
	}
	for (int j = 0; j < (int)E[S].size(); j++) {
		if (E[S][j]->getTo() == S) {
			continue;
		}
		V[ E[S][j]->getTo() ]->increaseExcess( E[S][j]->getG() );
		V[S]->increaseExcess(-E[S][j]->getG());
		E[S][j]->push(E[S][j]->getG());
	}
	V[S]->setLabel(V.size());
}

void FlowGraph::discharge(int u, int S, int T, queue<int>& q, vector<int>& used) {
	while (V[u]->getExcess() > 0) {
		if (V[u]->getHead() == (int)E[u].size()) {
			relabel(u);
			V[u]->initHead();
		}
		else {
			FlowEdge* e = E[u][ V[u]->getHead() ];
			if (e->getG() > 0 && V[u]->getLabel() == V[e->getTo()]->getLabel() + 1) {
				push(e);
				if (used[e->getTo()] == 0 && V[e->getTo()]->getExcess() > 0 && e->getTo() != S && e->getTo() != T) {
					q.push(e->getTo());
					used[e->getTo()] = 1;
				}
			}
			else
				V[u]->increaseHead();
		}
	}
	used[u] = 0;
}

long long FlowGraph::RelabelToFront(int S, int T) {
	initPreflow(S);
	vector<int> used(V.size());
	queue<int> q;
	for (int i = 0; i < (int)V.size(); i++) {
		V[i]->initHead();
		if (i != S && i != T && V[i]->getExcess() > 0) {
			q.push(i);
			used[i] = 1;
		}
	}
	while (q.size() > 0) {
		int cur = q.front();
		q.pop();

		discharge(cur, S, T, q, used);
	}

	return V[T]->getExcess();
}

void FlowGraph::printFlow() {
	vector<int> ans;
	ans.resize(EdgeCount);
	for (int i = 0; i < (int)E.size(); i++) {
		for (int j = 0; j < (int)E[i].size(); j++) {
			if (E[i][j]->getId() != 0) {
				int id = E[i][j]->getId();
				ans[id - 1] = E[i][j]->getFlow();
			}
		}
	}

	for (int i = 0; i < (int)ans.size(); i++) {
		cout << ans[i] << endl;
	}
}

int main() {
	FlowGraph G;
	G.get();
	cout << G.RelabelToFront(0, G.getSize() - 1) << endl;
	G.printFlow();
	return 0;
}
