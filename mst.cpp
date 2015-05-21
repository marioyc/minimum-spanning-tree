#include <cstdio>
#include <cstring>
#include <climits>
#include <algorithm>
#include <vector>

using namespace std;

#define MAXN 100
#define MAXM 10000

int n,m;

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

vector<int> L[MAXN],edge_id[MAXN],path,path_edges;
bool visited[MAXN];

bool dfs(int u, int v){
	path.push_back(u);

	if(u == v) return true;

	bool ret = false;

	if(!visited[u]){
		visited[u] = true;

		for(int i = 0;i < L[u].size();++i){
			int to = L[u][i];
			path_edges.push_back(edge_id[u][i]);

			if(dfs(to,v)){
				ret = true;
				break;
			}else{
				path_edges.pop_back();
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

// dfs to root the tree
int a[MAXN],in[MAXN],out[MAXN],cont;
int parent_lca[10][MAXN];
int height[MAXN];

void dfs3(int cur, int p = 0){
	a[cont] = cur;
	in[cur] = cont++;

	parent_lca[0][cur] = p;
	height[cur] = (p == 0? 0 : 1 + height[p]);

	for(int i = 1;i <= 9;++i)
		parent_lca[i][cur] = parent_lca[i - 1][ parent_lca[i - 1][cur] ];

	for(int i = (int)L[cur].size() - 1;i >= 0;--i){
		int to = L[cur][i];

		if(to != p)
			dfs3(to,cur);
	}

	out[cur] = cont;
}

// lowest common ancestor

int lca(int u, int v){
	if(height[u] < height[v])
		swap(u,v);

	for(int i = 9;i >= 0;--i)
		if((height[u] - height[v]) >> i & 1)
			u = parent_lca[i][u];
	
	if(u == v)
		return u;

	for(int i = 9;i >= 0;--i)
		if(parent_lca[i][u] != parent_lca[i][v]){
			u = parent_lca[i][u];
			v = parent_lca[i][v];
		}

	return parent_lca[0][u];
}

// dfs for the sensitivity analysis of tree edges

int k[MAXN],tree_sensitivity[MAXN];

void dfs4(int cur, int p){
	if(p != 0)
		tree_sensitivity[cur] = query(0,0,n - 1,in[cur],out[cur]);

	// update?

	for(int i = (int)L[cur].size() - 1;i >= 0;--i){
		int to = L[cur][i];

		if(to != p)
			dfs4(to,cur);
	}
}

int main(){
	n = m = 0;

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
			edge_id[ e[i].u ].push_back(i);

			L[ e[i].v ].push_back(e[i].u);
			edge_id[ e[i].v ].push_back(i);

			printf("%d - %d (%d)\n",e[i].u,e[i].v,e[i].w);
		}else in_tree[i] = false;
	}

	printf("Total weight = %d\n\n",total);

	// Find the path int T that connect two endpoint of an edge

	bool verified = true;

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

			// Verify MST

			for(int j = 0;j < path_edges.size();++j){
				if(e[i].w < e[ path_edges[j] ].w)
					verified = false;
			}
		}
	}

	printf("\n");

	printf("Verification: ");
	if(verified) printf("It is an MST\n\n");
	else printf("It is not an MST\n\n");

	// Find edges that can replace an edge in the MST

	int most_vital_edge = -1;
	int min_increase = -1;

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

			if(most_vital_edge == -1 || e[rep].w - e[i].w < min_increase){
				min_increase = e[rep].w - e[i].w;
				most_vital_edge = i;
			}

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
	printf("Most vital edge = (%d, %d)\n",e[most_vital_edge].u,e[most_vital_edge].v);
	printf("Increase = %d\n",min_increase);

	// Sensitivity analysis for tree edges

	cont = 0;
	dfs3(1,0);

	for(int i = 1;i <= n;++i)
		k[i] = INT_MAX;

	dfs4(1,0);

	return 0;
}