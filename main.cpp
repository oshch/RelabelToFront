#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <list>

using namespace std;

const int INF = 1e9;
const long long INFLL = (long long)1e18;

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

class Vertex {
private:
	int dist;
	list<Edge*>::iterator edges_head;
	long long excess;
	int label;
public:
	Vertex() : dist(INF), edges_head(NULL) {}

	int getDist();
	list<Edge*>::iterator getHead();
	void setDist(int x);
	void increaseHead();
	void initHead(list<Edge*>::iterator i);

	long long getExcess();
	int getLabel();
	void setExcess(long long x);
	void setLabel(int x);
	void increaseExcess(long long x);
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

class Graph {
protected:
	int EdgeCount;
	vector<list<Edge*> > E;
	vector<Vertex*> V;

public:
	Graph() {}
	Graph(vector<list<Edge*> >& FlowE, const vector<Vertex*>& FlowV, int minG);
	int getSize();
	void bfs(int S);
};

class FlowGraph : public Graph {
public:
	FlowGraph() {}

	int getMaxCap();
	void get();
	long long pushFlow(int S, long long flow, int minG, int T);
	long long Dinic(int S, int T);
	void printFlow();

	void push(FlowEdge* e);
	void relabel(int a);
	void initPreflow(int S);

	void discharge(int u, int S, int T, queue<int>& q, vector<int>& used);
	long long RelabelToFront(int S, int T);
};

int Vertex::getDist() {
	return dist;
}

list<Edge*>::iterator Vertex::getHead() {
	return edges_head;
}

void Vertex::setDist(int x) {
	dist = x;
}

void Vertex::increaseHead() {
	edges_head++;
}

void Vertex::initHead(list<Edge*>::iterator i) {
	edges_head = i;
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

void Graph::bfs(int S) {
	int n = E.size();
	for (int i = 0; i < n; i++) {
		V[i]->setDist(INF);
	}
	V[S]->setDist(0);
	queue<int> q;
	q.push(S);

	while (q.size() > 0) {
		int cur = q.front();
		q.pop();

		for (list<Edge*>::iterator i = E[cur].begin(); i != E[cur].end(); i++) {
			int y = (*i)->getTo();
			if (V[y]->getDist() == INF) {
				V[y]->setDist(V[cur]->getDist() + 1);
				q.push(y);
			}
		}
	}
}

Graph::Graph(vector<list<Edge*> >& FlowE, const vector<Vertex*>& FlowV, int minG) {
	E.resize(FlowE.size());
	EdgeCount = 0;
	for (int i = 0; i < (int)FlowE.size(); i++) {
		E[i].resize(0);
		for (list<Edge*>::iterator j = FlowE[i].begin(); j != FlowE[i].end(); j++) {
			if (((FlowEdge*) *j)->getG() >= minG) {
				int y = (*j)->getTo();
				E[i].push_back(new Edge(i, y, (*j)->getId()));
				EdgeCount++;
			}
		}
	}

	V = FlowV;
}

int FlowGraph::getMaxCap() {
	int ans = 0;
	for (int i = 0; i < (int)E.size(); i++) {
		for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); j++) {
			ans = max(ans, ((FlowEdge*) *j)->getFlow() + ((FlowEdge*) *j)->getG());
		}
	}
	return ans;
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
		E[b].push_back(new FlowEdge(b, a, 0, 0, 0, (FlowEdge*) E[a].back()));
		((FlowEdge*) E[a].back())->setRev((FlowEdge*) E[b].back());
	}
}

long long FlowGraph::pushFlow(int x, long long flow, int minG, int T) {
	if (x == T) {
		return flow;
	}
	long long oldflow = flow;

	for (list<Edge*>::iterator i = V[x]->getHead(); i != E[x].end(); i++) {
		int y = (*i)->getTo();
		if (V[y]->getDist() == V[x]->getDist() + 1 && ((FlowEdge*) *i)->getG() >= minG && flow > 0) {
			long long cur = pushFlow(y, min(flow, (long long)((FlowEdge*) *i)->getG()), minG, T);
			flow -= cur;
			if (flow != 0) {
				V[x]->increaseHead();
			}
			((FlowEdge*) *i)->push(cur);
		}
	}

	return (oldflow - flow);
}

long long FlowGraph::Dinic(int S, int T) {
	long long flow = 0;
	for (int minG = this->getMaxCap(); minG > 0; minG /= 2) {
		while (true) {
			Graph G(E, V, minG);
			G.bfs(S);
			for (int i = 0; i < (int)V.size(); i++) {
				V[i]->initHead(E[i].begin());
			}
			long long cur = pushFlow(S, INFLL, minG, T);
			if (cur == 0) {
				break;
			}
			flow += cur;
		}
	}
	return flow;
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
	for (list<Edge*>::iterator i = E[a].begin(); i != E[a].end(); i++) {
		if (((FlowEdge*) *i)->getG() > 0) {
			ans = min(ans, V[ (*i)->getTo() ]->getLabel() + 1);
		}
	}
	V[a]->setLabel(ans);
}

void FlowGraph::initPreflow(int S) {
	for (int i = 0; i < (int)V.size(); i++) {
		V[i]->setExcess(0);
		V[i]->setLabel(0);
	}
	for (list<Edge*>::iterator i = E[S].begin(); i != E[S].end(); i++) {
		if ((*i)->getTo() == S) {
			continue;
		}
		V[ (*i)->getTo() ]->increaseExcess( ((FlowEdge*) *i)->getG() );
		V[S]->increaseExcess(-((FlowEdge*) *i)->getG());
		((FlowEdge*) *i)->push(((FlowEdge*) *i)->getG());
	}
	V[S]->setLabel(V.size());
}

void FlowGraph::discharge(int u, int S, int T, queue<int>& q, vector<int>& used) {
	while (V[u]->getExcess() > 0) {
		if (V[u]->getHead() == E[u].end()) {
			relabel(u);
			V[u]->initHead(E[u].begin());
		}
		else {
			FlowEdge* e = (FlowEdge*) *V[u]->getHead();
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
		V[i]->initHead(E[i].begin());
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
		for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); j++) {
			if ((*j)->getId() != 0) {
				int id = (*j)->getId();
				ans[id - 1] = ((FlowEdge*) *j)->getFlow();
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
