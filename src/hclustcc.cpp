#include <Rcpp.h>
using namespace Rcpp;

// iPair ==> Integer Pair
typedef std::pair<int, float> iPair;

struct node
{
  int id;
  int size;
  NumericVector x;
  int father;
  float height;
  std::map<int,float,std::greater<float>> neibs;
};

float dist(NumericMatrix::Row x1,NumericMatrix::Row x2){
  return sum((x1-x2)*(x1-x2));
}


float dist(NumericVector x1,NumericVector x2){
  return sum((x1-x2)*(x1-x2));
}

float dist(node n1,node n2,String method){
  float d=0;
  if(method=="ward"){
    d=static_cast< float >(n1.size*n2.size)/static_cast< float >(n1.size+n2.size)*dist(n1.x,n2.x);
  }else if(method=="average" || method=="median"){
    d=dist(n1.x,n2.x);
  }else{
    stop("Unknown method");
  }
  //return 
  return d;
}


void print_pq(std::multimap<float,std::pair<int, int>,std::less<float>> priority_queue){
  for(auto it=priority_queue.begin();it!=priority_queue.end();it++){
    std::pair<int,int> edge = it->second;
    Rcout << edge.first << " - " << edge.second << ":" << it->first << std::endl;
  }
}

// [[Rcpp::export]]
List hclustcc_cpp(List nb,NumericMatrix X,String method) {

  int V = X.nrow();

  
  // data-structure creation
  std::vector<node> graph(2*V-1);
  std::multimap<float,std::pair<int, int>,std::less<float>> priority_queue;
  for(int i=0; i<nb.length(); ++i){
    if(nb[i]!=R_NilValue) {
      NumericVector nbi = as<NumericVector>(nb[i]);
      NumericMatrix::Row x1 = X(i,_);
      node cnode;
      cnode.x = x1;
      cnode.id = i;
      cnode.size=1;
      cnode.height=0;
      for(int n=0; n<nbi.length(); ++n){
        int j = nbi[n];
        if(i!=j){
          NumericMatrix::Row x2 = X(j,_);
          float d = dist(x1,x2);
          cnode.neibs.insert(std::make_pair(j,d));
          if(i<j){
            priority_queue.insert(std::make_pair(d,std::make_pair(i,j)));
          }
        }
      }
      graph[i]=cnode;
    }
  }
  
  // Lets Merge !
  float H =0;
  NumericMatrix merge(V-1,2);
  NumericVector height(V-1);
  
  for(int imerge=0;imerge<(V-1);imerge++){
    

    int node_id = V+imerge;
    auto best_merge = priority_queue.begin();
    // Rcout << "###############"  << std::endl;
    // print_pq(priority_queue);
    if(best_merge==priority_queue.end()){
      // deal with isolated regions (no more merge possible and node)
      // Rcout << "isolated regions" << std::endl;
      for(int nb_miss_merge=node_id;nb_miss_merge<(V-1);nb_miss_merge++){
        merge(nb_miss_merge,1)=0;
        merge(nb_miss_merge,2)=0;
      }
      break;
    }
    std::pair<int,int> edge = best_merge->second; 
    int g = std::get<0>(edge);
    int h = std::get<1>(edge);
    
    H = H + best_merge->first;
    height[imerge]=sqrt(best_merge->first);
    
    // strore merge move in hclust format with 1 based indices
    if(g<V){
      merge(imerge,0)=-(g+1);
    }else{
      merge(imerge,0)=g-V+1;
    }

    if(h<V){
      merge(imerge,1)=-(h+1);
    }else{
      merge(imerge,1)=h-V+1;
    }
    
    // create a new node
    node new_node;
    
    // update the stored position
    if(method=="ward" || method=="average"){
        new_node.x  = (graph[g].size*graph[g].x + graph[h].size*graph[h].x)/(graph[g].size+graph[h].size);      
    }else if(method=="median"){
        new_node.x  = (graph[g].x + graph[h].x)/2;   
    }else{
        stop("Unknown method");
    }

    new_node.id = node_id;
    new_node.size=graph[g].size+graph[h].size;
    new_node.height = H;
    new_node.father=0;
    node node_g = graph[g];
    node node_h = graph[h];
    // stockage de la fusion
    graph[g].father = node_id;
    graph[h].father = node_id;
    

    for(auto nei_g = node_g.neibs.begin();nei_g!=node_g.neibs.end();nei_g++){
      
      int i = g;
      int j = nei_g->first;
      float v = nei_g->second;
      
      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(v);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          s=priority_queue.erase(s);
          if(s==search.second){
            break;
          }
        }
      }
      
      // old link deletion in graph
      graph[j].neibs.erase(i);
      // new links in graph
      if(j!=h){
        // distance calculation
        float d = dist(new_node,graph[j],method);
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
      
    }
    
    
    for(auto nei_h = node_h.neibs.begin();nei_h!=node_h.neibs.end();nei_h++){
      
      int i = h;
      int j = nei_h->first;
      float v = nei_h->second;
      
      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(v);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          s=priority_queue.erase(s);
          if(s==search.second){
            break;
          }
        }
      }
      // old link deletion in graph
      graph[j].neibs.erase(i);
      
      // new links in graphs
      if(j!=g){
        float d = dist(new_node,graph[j],method);
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
    }
    
    
    // add the newly created node
    graph[node_id]=new_node;
    
    // add the new possible merges in the priority queue
    for(auto nei = new_node.neibs.begin();nei!=new_node.neibs.end();nei++){
      float d = nei->second;
      float j = nei->first;
      priority_queue.insert(std::make_pair(d,std::make_pair(j,node_id)));
    }

    
  }
  
  NumericVector tree(2*V-1);
  NumericVector tree_height(2*V-1);
  NumericVector ids(2*V-1);
  for(int i=0;i<(2*V-1);i++){
    node cnode = graph[i];
    tree(i)=cnode.father+1;
    tree_height(i)=cnode.height;
    ids(i)=cnode.id+1;
  }
  
  // 
  return List::create(Named("tree",tree),Named("H",tree_height),Named("id",ids),Named("merge",merge),Named("height",height));
  
}


