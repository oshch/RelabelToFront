#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <stack>
#include <list>

using namespace std;

const int INF = 1e9;
const long long INFLL = (long long)1e16;

class Edge {
protected:
	int from, to;
	int id;
	bool used;
public:
	Edge() {}
	Edge(int _from, int _to, int _id) : from(_from), to(_to), id(_id) {}
	void setUsed(bool u);
	bool getUsed();
	int getTo();
	int getId();
	int getFrom();
};

class Vertex {
protected:
	int dist;
	list<Edge*>::iterator edges_head;
	long long excess;
	int label, id;
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
	void setId(int x);
	int getId();
};

class FlowEdge : public Edge {
protected:
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
	int getCapacity();
	FlowEdge* getRev();
	void setFlow(int f);
};

class MKMVertex : public Vertex {
private:
    long long inMaxFlow, outMaxFlow;
    bool deleted;
public:
    MKMVertex(){reset();}
    MKMVertex(Vertex* v) {
        id = v->getId();
        reset();
    }
    void addInFlow(int df);
    void addOutFlow(int df);
    long long getPotential();
    void reset();
    void setDeleted(bool b);
    bool isDeleted();
};

class MKMEdge : public FlowEdge {
private:
    FlowEdge* parent;
    list<Edge*>::iterator inIter, outIter;
public:
    MKMEdge(FlowEdge* edge) {
        parent = edge;
        from = edge->getFrom();
        to = edge->getTo();
        id = edge->getId();
        flow = edge->getFlow();
        capacity = edge->getCapacity();
    }
    long long pushFlow(MKMVertex* from, MKMVertex* to);
    FlowEdge* getParent();
    void setInIter(list<Edge*>::iterator i);
    void setOutIter(list<Edge*>::iterator i);
    list<Edge*>::iterator getInIter();
    list<Edge*>::iterator getOutIter();
};

class Graph {
friend class MKMFinder;
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

class MKMFinder : public Graph {
private:
    vector<list<Edge*> > inE;
    vector<long long> potential;
    FlowGraph* parent;
    int S, T;

    void removeEdge(MKMEdge* edge, int from, int to);
	void pushFlow(vector<list<Edge*> >& edges, int S);
public:
	MKMFinder() {}
	MKMFinder(FlowGraph* graph, int S, int T) : S(S), T(T) {
        for (int i = 0; i < graph->getSize(); ++i) {
            V.push_back(new MKMVertex(graph->V[i]));
        }
        E.assign(graph->getSize(), list<Edge*>());
        inE.assign(graph->getSize(), list<Edge*>());

        int n = graph->E.size();
        for (int i = 0; i < n; ++i) {
            V[i]->setDist(INF);
        }
        V[S]->setDist(0);
        queue<int> q;
        q.push(S);

        while (!q.empty()) {
            int cur = q.front();
            q.pop();

            for (list<Edge*>::iterator i = graph->E[cur].begin(); i != graph->E[cur].end(); ++i) {
                if (((FlowEdge*) (*i))->getCapacity() - ((FlowEdge*) (*i))->getFlow() <= 0)
                    continue;
                int y = (*i)->getTo();
                if (V[y]->getDist() == INF) {
                    V[y]->setDist(V[cur]->getDist() + 1);
                    (*i)->setUsed(true);
                    q.push(y);
                }
            }
        }

        for (int i = 0; i < graph->getSize(); ++i) {
            for(list<Edge*>::iterator j = graph->E[i].begin(); j != graph->E[i].end(); ++j) {
                if (V[(*j)->getFrom()]->getDist() != V[(*j)->getTo()]->getDist() - 1)
                    continue;
                E[i].push_back((Edge*) new MKMEdge((FlowEdge*) *j));
                inE[(*j)->getTo()].push_back(E[i].back());
                ((MKMEdge*) E[i].back())->setInIter(--(inE[(*j)->getTo()].end()));
                ((MKMEdge*) E[i].back())->setOutIter((--E[i].end()));
            }
        }
	}
	~MKMFinder() {
        for (int i = 0; i < getSize(); ++i) {
            delete V[i];
        }
        for (int i = 0; i < getSize(); ++i) {
            for(list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
                delete (*j);
            }
        }
	}

    long long solve();
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

void Vertex::setId(int x) {
    id = x;
}

int Vertex::getId() {
    return id;
}

void Edge::setUsed(bool u) {
    used = u;
}
bool Edge::getUsed() {
    return used;
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

int FlowEdge::getCapacity() {
    return capacity;
}
FlowEdge* FlowEdge::getRev() {
    return rev;
}
void FlowEdge::setFlow(int f) {
    flow = f;
}

void MKMVertex::addInFlow(int df) {
    inMaxFlow += df;
}
void MKMVertex::addOutFlow(int df) {
    outMaxFlow += df;
}
long long MKMVertex::getPotential() {
    return min(inMaxFlow, outMaxFlow);
}
void MKMVertex::reset() {
    inMaxFlow = outMaxFlow = excess = 0;
    deleted = false;
}
void MKMVertex::setDeleted(bool b) {
    deleted = b;
}
bool MKMVertex::isDeleted() {
    return deleted;
}

long long MKMEdge::pushFlow(MKMVertex* from, MKMVertex* to) {
    long long f = min(from->getExcess(), min(to->getPotential(), 0ll + capacity - flow));
    if (from->getId() == getFrom()) {
        from->addOutFlow(-f);
        to->addInFlow(-f);
    } else {
        to->addOutFlow(-f);
        from->addInFlow(-f);
    }
    flow += f;
    to->increaseExcess(f);
    from->increaseExcess(-f);
    return f;
}
FlowEdge* MKMEdge::getParent() {
    return parent;
}
void MKMEdge::setInIter(list<Edge*>::iterator i) {
    inIter = i;
}
void MKMEdge::setOutIter(list<Edge*>::iterator i) {
    outIter = i;
}
list<Edge*>::iterator MKMEdge::getInIter() {
    return inIter;
}
list<Edge*>::iterator MKMEdge::getOutIter() {
    return outIter;
}

int Graph::getSize() {
	return E.size();
}

void Graph::bfs(int S) {
	int n = E.size();
	for (int i = 0; i < n; ++i) {
		V[i]->setDist(INF);
	}
	V[S]->setDist(0);
	queue<int> q;
	q.push(S);

	while (q.size() > 0) {
		int cur = q.front();
		q.pop();

		for (list<Edge*>::iterator i = E[cur].begin(); i != E[cur].end(); ++i) {
            (*i)->setUsed(false);
			int y = (*i)->getTo();
			if (V[y]->getDist() == INF) {
				V[y]->setDist(V[cur]->getDist() + 1);
				(*i)->setUsed(true);
				q.push(y);
			}
		}
	}
}

Graph::Graph(vector<list<Edge*> >& FlowE, const vector<Vertex*>& FlowV, int minG) {
	E.resize(FlowE.size());
	EdgeCount = 0;
	for (int i = 0; i < (int)FlowE.size(); ++i) {
		E[i].resize(0);
		for (list<Edge*>::iterator j = FlowE[i].begin(); j != FlowE[i].end(); ++j) {
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
	for (int i = 0; i < (int)E.size(); ++i) {
		for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
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
	for (int i = 0; i < n; ++i) {
		V[i] = new Vertex();
		V[i]->setId(i);
	}

	for (int i = 0; i < m; ++i) {
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

	for (list<Edge*>::iterator i = V[x]->getHead(); i != E[x].end(); ++i) {
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
    while (true) {
        MKMFinder net(this, S, T);
        long long cur = net.solve();
        if (cur == 0) {
            break;
        }
        flow += cur;
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
	for (list<Edge*>::iterator i = E[a].begin(); i != E[a].end(); ++i) {
		if (((FlowEdge*) *i)->getG() > 0) {
			ans = min(ans, V[ (*i)->getTo() ]->getLabel() + 1);
		}
	}
	V[a]->setLabel(ans);
}

void FlowGraph::initPreflow(int S) {
	for (int i = 0; i < (int)V.size(); ++i) {
		V[i]->setExcess(0);
		V[i]->setLabel(0);
	}
	for (list<Edge*>::iterator i = E[S].begin(); i != E[S].end(); ++i) {
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
	for (int i = 0; i < (int)V.size(); ++i) {
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
	for (int i = 0; i < (int)E.size(); ++i) {
		for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
			if ((*j)->getId() != 0) {
				int id = (*j)->getId();
				ans[id - 1] = ((FlowEdge*) *j)->getFlow();
			}
		}
	}

	for (int i = 0; i < (int)ans.size(); ++i) {
		cout << ans[i] << endl;
	}
}

void MKMFinder::removeEdge(MKMEdge* edge, int from, int to) {
    long long df = edge->getCapacity() - edge->getFlow();
    inE[to].erase(edge->getInIter());
    E[from].erase(edge->getOutIter());
    ((MKMVertex*) V[from])->addOutFlow(-df);
    ((MKMVertex*) V[to])->addInFlow(-df);
    edge->getParent()->setFlow(edge->getFlow());
    edge->getParent()->getRev()->setFlow(-edge->getFlow());
    delete edge;
}

void MKMFinder::pushFlow(vector<list<Edge*> >& edges, int S) {
    queue<int> q;
    q.push(S);
    while(!q.empty()) {
        int currNum = q.front();
        q.pop();
        MKMVertex* current = ((MKMVertex*) V[currNum]);
        for (list<Edge*>::iterator i = edges[currNum].begin(); i != edges[currNum].end(); ++i) {
            if (current->getExcess() == 0)
                break;
            int to = (*i)->getTo();
            if (to == currNum)
                to = (*i)->getFrom();
            ((MKMEdge*) (*i))->pushFlow((MKMVertex*) V[currNum], (MKMVertex*) V[to]);
            q.push(to);
        }
    }
}

long long MKMFinder::solve() {
    for (int i = 0; i < getSize(); ++i) {
        ((MKMVertex *) V[i])->reset();
    }
    for (int i = 0; i < getSize(); ++i) {
        for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
            ((MKMVertex *) V[(*j)->getFrom()])->addOutFlow(((MKMEdge*) (*j))->getCapacity() - ((MKMEdge*) (*j))->getFlow());
            ((MKMVertex *) V[(*j)->getTo()])->addInFlow(((MKMEdge*) (*j))->getCapacity() - ((MKMEdge*) (*j))->getFlow());
        }
    }
    ((MKMVertex *) V[S])->addInFlow(INFLL);
    ((MKMVertex *) V[T])->addOutFlow(INFLL);
    stack<int> vertices;
    long long totalFlow = 0;
    while (true) {
        for (int i = 0; i < getSize(); ++i) {
            if (((MKMVertex*) V[i])->getPotential() == 0 && !((MKMVertex*) V[i])->isDeleted()) {
                vertices.push(i);
            }
        }
        while (!vertices.empty()) {
            int id = vertices.top();
            vertices.pop();
            ((MKMVertex*) V[id])->setDeleted(true);
            while (!inE[id].empty()) {
                MKMEdge* e = (MKMEdge*) *inE[id].begin();
                int from = e->getFrom();
                int to = e->getTo();
                removeEdge(e, from, to);
                if (((MKMVertex*) V[from])->getPotential() == 0 && !((MKMVertex*) V[from])->isDeleted())
                    vertices.push(from);
                if (((MKMVertex*) V[to])->getPotential() == 0 && !((MKMVertex*) V[to])->isDeleted())
                    vertices.push(to);
            }
            while (!E[id].empty()) {
                MKMEdge* e = (MKMEdge*) *E[id].begin();
                int from = e->getFrom();
                int to = e->getTo();
                removeEdge(e, from, to);
                if (((MKMVertex*) V[from])->getPotential() == 0 && !((MKMVertex*) V[from])->isDeleted())
                    vertices.push(from);
                if (((MKMVertex*) V[to])->getPotential() == 0 && !((MKMVertex*) V[to])->isDeleted())
                    vertices.push(to);
            }

        }
        long long minPotential = INFLL;
        MKMVertex* minVertex;
        for (int i = 0; i < getSize(); ++i) {
            if (!((MKMVertex*) V[i])->isDeleted() && minPotential > ((MKMVertex*) V[i])->getPotential()) {
                minVertex = ((MKMVertex*) V[i]);
                minPotential = minVertex->getPotential();
            }
        }
        if (minPotential == INFLL)
            break;
        minVertex->increaseExcess(minPotential);
        //print();
        pushFlow(E, minVertex->getId());
        //print();
        minVertex->increaseExcess(minPotential);
        pushFlow(inE, minVertex->getId());
        //print();
        totalFlow += minPotential;
    }
    for (int i = 0; i < getSize(); ++i) {
        for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
            ((MKMEdge*) (*j))->getParent()->setFlow(((MKMEdge*) (*j))->getFlow());
            ((MKMEdge*) (*j))->getParent()->getRev()->setFlow(-((MKMEdge*) (*j))->getFlow());
        }
    }
    return totalFlow;
}


int main() {
    //freopen("input.txt", "rt", stdin);
	FlowGraph G;
	G.get();
	cout << G.Dinic(0, G.getSize() - 1) << endl;
	G.printFlow();
	return 0;
}
