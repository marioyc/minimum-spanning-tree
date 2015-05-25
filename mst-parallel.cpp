#include <cstdio>
#include <cstring>
#include <climits>
#include <algorithm>
#include <vector>

#include <mpi.h>

using namespace std;

#define MAXN 100
#define MAXM 10000

int p,rank;
int n,m;

struct edge{
	int u,v,w;
	int id;

	bool operator < (edge e) const{
		return w < e.w;
	}
}e[MAXM],e_local[MAXM];

bool in_tree[MAXM];

// Union - Find

int parent[MAXN + 1],uf_rank[MAXN + 1];
int Find(int x){
	if(parent[x] == x) return x;
	parent[x] = Find(parent[x]);
	return parent[x];
}

void Union(int x, int y){
	x = Find(x);
	y = Find(y);

	if(uf_rank[x] < uf_rank[y]){
		parent[x] = y;
	}else{
		parent[y] = x;

		if(uf_rank[x] == uf_rank[y])
			++uf_rank[x];
	}
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

int main(int argc, char *argv[]){
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(rank == 0){
		n = m = 0;

		while(scanf("%d,%d,%d",&e[m].u,&e[m].v,&e[m].w) == 3){
			n = max(n,max(e[m].u,e[m].v));
			e[m].id = m;
			++m;
		}
	}

	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD);

	int blocklengths[] = {1,1,1,1};
	MPI_Aint displacements[] = {offsetof(edge,u), offsetof(edge,v), offsetof(edge,w), offsetof(edge,id)};
	MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype MPI_EDGE;

	MPI_Type_create_struct(4,blocklengths,displacements,types,&MPI_EDGE);
	MPI_Type_commit(&MPI_EDGE);

	MPI_Bcast(e,m,MPI_EDGE,0,MPI_COMM_WORLD);

	// Find minimum spanning tree using Kruskal's algorithm

	for(int i = 1;i <= n;++i){
		parent[i] = i;
		uf_rank[i] = 0;
	}

	int start = (m + p - 1) / p * rank,end = min(m,(m + p - 1) / p * (rank + 1));
	int m_local = 0;

	for(int i = start;i < end;++i)
		e_local[i - start] = e[i];

	sort(e_local,e_local + (end - start));

	for(int i = 0;i < end - start;++i){
		if(Find(e_local[i].u) != Find(e_local[i].v)){
			Union(e_local[i].u,e_local[i].v);
			e_local[m_local] = e_local[i];
			m_local++;
		}
	}

	int curP = p;

	while(curP > 1){
		if(rank < curP){
			//printf("%d %d\n",rank,curP);
			int m1 = 0,m2 = 0;
			MPI_Status status;

			//printf("%d %d send\n",rank,curP);
			MPI_Send(&m_local,1,MPI_INT,rank / 2,0,MPI_COMM_WORLD);
			//printf("%d %d send\n",rank,curP);
			MPI_Send(e_local,m_local,MPI_EDGE,rank / 2,1,MPI_COMM_WORLD);
			//printf("%d %d send\n",rank,curP);

			if(2 * rank < curP){
				//printf("%d %d receive\n",rank,curP);
				MPI_Recv(&m1,1,MPI_INT,2 * rank,0,MPI_COMM_WORLD,&status);
				//printf("%d %d receive\n",rank,curP);
				MPI_Recv(e_local,m1,MPI_EDGE,2 * rank,1,MPI_COMM_WORLD,&status);
				//printf("%d %d receive\n",rank,curP);


				if(2 * rank + 1 < curP){
					//printf("%d %d receive\n",rank,curP);
					MPI_Recv(&m2,1,MPI_INT,2 * rank + 1,0,MPI_COMM_WORLD,&status);
					//printf("%d %d receive\n",rank,curP);
					MPI_Recv(e_local + m1,m2,MPI_EDGE,2 * rank + 1,1,MPI_COMM_WORLD,&status);
					//printf("%d %d receive\n",rank,curP);
				}

				int m_aux = m1 + m2;
				sort(e_local,e_local + m_aux);

				for(int i = 1;i <= n;++i){
					parent[i] = i;
					uf_rank[i] = 0;
				}

				m_local = 0;

				for(int i = 0;i < m_aux;++i){
					if(Find(e_local[i].u) != Find(e_local[i].v)){
						Union(e_local[i].u,e_local[i].v);
						e_local[m_local] = e_local[i];
						m_local++;
					}
				}
			}
		}

		curP = (curP + 1) / 2;
	}

	if(rank == 0){
		int total = 0;

		for(int i = 0;i < m_local;++i){
			L[ e_local[i].u ].push_back(e_local[i].v);
			edge_id[ e_local[i].u ].push_back(e_local[i].id);

			L[ e_local[i].v ].push_back(e_local[i].u);
			edge_id[ e_local[i].v ].push_back(e_local[i].id);

			in_tree[ e_local[i].id ] = true;

			printf("%d - %d (%d)\n",e_local[i].u,e_local[i].v,e_local[i].w);
			total += e_local[i].w;
		}

		printf("Total weight = %d\n\n",total);
	}

	for(int i = 1;i <= n;++i){
		int list_size;

		if(rank == 0)
			list_size = L[i].size();

		MPI_Bcast(&list_size,1,MPI_INT,0,MPI_COMM_WORLD);

		L[i].resize(list_size);
		MPI_Bcast(&L[i].front(),list_size,MPI_INT,0,MPI_COMM_WORLD);

		edge_id[i].resize(list_size);
		MPI_Bcast(&edge_id[i].front(),list_size,MPI_INT,0,MPI_COMM_WORLD);
	}

	MPI_Bcast(in_tree,m,MPI_C_BOOL,0,MPI_COMM_WORLD);

	// Find the path int T that connect two endpoint of an edge

	bool verified_local = true;

	for(int i = (m + p - 1) / p * rank;i < min(m,(m + p - 1) / p * (rank + 1));++i){
		if(!in_tree[i]){
			memset(visited,false,sizeof visited);
			path.clear();
			dfs(e[i].u,e[i].v);

			printf("process %d, (%d, %d) : ",rank,e[i].u,e[i].v);

			for(int j = 0;j < path.size();++j){
				if(j > 0) printf(" -> ");
				printf("%d",path[j]);
			}

			printf("\n");

			// Verify MST

			for(int j = 0;j < path_edges.size();++j){
				if(e[i].w < e[ path_edges[j] ].w)
					verified_local = false;
			}
		}
	}

	printf("\n");

	bool verified_global;
	MPI_Reduce(&verified_local,&verified_global,1,MPI_C_BOOL,MPI_LAND,0,MPI_COMM_WORLD);

	if(rank == 0){
		printf("Verification: ");
		if(verified_global) printf("It is an MST\n\n");
		else printf("It is not an MST\n\n");
	}

	// Find edges that can replace an edge in the MST

	int most_vital_edge_local[p];
	int max_increase_local[p];
	most_vital_edge_local[rank] = -1;
	max_increase_local[rank] = -1;

	for(int i = (m + p - 1) / p * rank;i < min(m,(m + p - 1) / p * (rank + 1));++i){
		if(in_tree[i]){
			memset(visited,false,sizeof visited);
			dfs2(e[i].u,e[i].v,e[i].u);
			dfs2(e[i].v,e[i].u,e[i].v);

			printf("process %d, (%d, %d) :",rank,e[i].u,e[i].v);

			int rep = -1;

			for(int j = 0;j < m;++j)
				if(j != i && ((mark[ e[j].u ] == e[i].u && mark[ e[j].v ] == e[i].v) || (mark[ e[j].v ] == e[i].u && mark[ e[j].u ] == e[i].v))){
					printf(" (%d, %d)",e[j].u,e[j].v);

					if(rep == -1 || e[j].w < e[rep].w)
						rep = j;
				}

			printf(", replacement : (%d, %d)\n",e[rep].u,e[rep].v);

			if(e[rep].w - e[i].w > max_increase_local[rank]){
				max_increase_local[rank] = e[rep].w - e[i].w;
				most_vital_edge_local[rank] = i;
			}
		}
	}

	if(rank != 0){
		MPI_Send(&max_increase_local[rank],1,MPI_INT,0,2,MPI_COMM_WORLD);
		MPI_Send(&most_vital_edge_local[rank],1,MPI_INT,0,3,MPI_COMM_WORLD);
	}else{
		int most_vital_edge = most_vital_edge_local[0],max_increase = max_increase_local[0];
		MPI_Status status;

		for(int i = 1;i < p;++i){
			MPI_Recv(&max_increase_local[i],1,MPI_INT,i,2,MPI_COMM_WORLD,&status);
			MPI_Recv(&most_vital_edge_local[i],1,MPI_INT,i,3,MPI_COMM_WORLD,&status);

			if(max_increase_local[i] > max_increase){
				max_increase = max_increase_local[i];
				most_vital_edge = most_vital_edge_local[i];
			}
		}

		printf("\n");
		printf("Most vital edge = (%d, %d)\n",e[most_vital_edge].u,e[most_vital_edge].v);
		printf("Increase = %d\n",max_increase);
	}

	MPI_Finalize();

	return 0;
}