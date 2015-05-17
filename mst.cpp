#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>

using namespace std;

#define MAXN 100
#define MAXM 10000

struct edge{
	int u,v,w;

	bool operator < (edge e) const{
		return w < e.w;
	}
}e[MAXM];

bool in_tree[MAXM];

// Union - Find

int parent[MAXN + 1];

int Find(int x){
	if(parent[x] == x) return x;
	parent[x] = Find(parent[x]);
	return parent[x];
}

void Union(int x, int y){
	x = Find(x);
	y = Find(y);
	parent[x] = y;
}

// dfs to find a path between two nodes in the MST

vector<int> L[MAXN],path;
bool visited[MAXN];

bool dfs(int u, int v){
	path.push_back(u);

	if(u == v) return true;

	bool ret = false;

	if(!visited[u]){
		visited[u] = true;

		for(int i = 0;i < L[u].size();++i){
			int to = L[u][i];

			if(dfs(to,v)){
				ret = true;
				break;
			}
		}
	}

	if(!ret)
		path.pop_back();

	return ret;
}

// dfs to mark edges in a subtree after erasing an edge

int mark[MAXN];

void dfs2(int u, int v, int cur){
	if(visited[cur]) return;
	visited[cur] = true;

	mark[cur] = u;

	for(int i = 0;i < L[cur].size();++i){
		int to = L[cur][i];

		if(to != u && to != v)
			dfs2(u,v,to);
	}
}

int main(){
	int n = 0,m = 0;

	while(scanf("%d,%d,%d",&e[m].u,&e[m].v,&e[m].w) == 3){
		n = max(n,max(e[m].u,e[m].v));
		++m;
	}

	// Find minimum spanning tree using Kruskal's algorithm

	for(int i = 1;i <= n;++i)
		parent[i] = i;

	sort(e,e + m);

	int total = 0;

	for(int i = 0;i < m;++i){
		if(Find(e[i].u) != Find(e[i].v)){
			Union(e[i].u,e[i].v);
			total += e[i].w;
			in_tree[i] = true;

			L[ e[i].u ].push_back(e[i].v);
			L[ e[i].v ].push_back(e[i].u);

			printf("%d - %d (%d)\n",e[i].u,e[i].v,e[i].w);
		}else in_tree[i] = false;
	}

	printf("Total weight = %d\n\n",total);

	// Find the path int T that connect two endpoint of an edge

	for(int i = 0;i < m;++i){
		if(!in_tree[i]){
			memset(visited,false,sizeof visited);
			path.clear();
			dfs(e[i].u,e[i].v);

			printf("(%d, %d) : ",e[i].u,e[i].v);

			for(int j = 0;j < path.size();++j){
				if(j > 0) printf(" -> ");
				printf("%d",path[j]);
			}

			printf("\n");
		}
	}

	printf("\n");

	// Find edges that can replace an edge in the MST

	for(int i = 0;i < m;++i){
		if(in_tree[i]){
			memset(visited,false,sizeof visited);
			dfs2(e[i].u,e[i].v,e[i].u);
			dfs2(e[i].v,e[i].u,e[i].v);

			printf("(%d, %d) :",e[i].u,e[i].v);

			int rep = -1;

			for(int j = 0;j < m;++j)
				if(j != i && ((mark[ e[j].u ] == e[i].u && mark[ e[j].v ] == e[i].v) || (mark[ e[j].v ] == e[i].u && mark[ e[j].u ] == e[i].v))){
					printf(" (%d, %d)",e[j].u,e[j].v);

					if(rep == -1 || e[j].w < e[rep].w)
						rep = j;
				}

			printf(", replacement : (%d, %d)\n",e[rep].u,e[rep].v);

			/*
			//test
			
			for(int j = 1;j <= n;++j)
				parent[j] = j;

			for(int j = 0;j < m;++j){
				if(j != i && Find(e[j].u) != Find(e[j].v)){
					Union(e[j].u,e[j].v);
					if(!in_tree[j])
						printf("(%d, %d)\n",e[j].u,e[j].v);
				}
			}*/
		}
	}

	printf("\n");

	return 0;
}