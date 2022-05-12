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
  float d;
  float pi;
  float L;
  float Lt;
  float r;
  std::map<int,float,std::greater<float>> neibs;
};


float log_dirichlet_multinom(NumericVector x){
  float beta = 1;
  int d = x.length();
  float icl_emiss = lgamma(d*beta)+sum(lgamma(x+beta))-d*lgamma(beta)-lgamma(sum(x+beta));
  return icl_emiss;
}

float ladd(float a, float b){
  float res = 0;
  if(a>b){
    res = a+log(1+exp(b-a));
  }else{
    res = b+log(1+exp(a-b));
  }
  return res;
}

float dist(NumericVector x1,NumericVector x2){
  return sum((x1-x2)*(x1-x2));
}

float distchi2(NumericVector x1,NumericVector x2,NumericVector w){
  if(sum(x1)==0||sum(x2)==0){
    return 0;
  }
  NumericVector x1n = x1 / sum(x1);
  NumericVector x2n = x2 / sum(x2);
  return sum(w*(x1n-x2n)*(x1n-x2n));
}

float dist(node n1,node n2,NumericVector w,String method){
  float d=0;
  if(method=="ward"){
    d=static_cast< float >(n1.size*n2.size)/static_cast< float >(n1.size+n2.size)*dist(n1.x,n2.x);
  }else if(method=="centroid" || method=="median"){
    d=dist(n1.x,n2.x);
  }else if(method=="chi2"){
    d=distchi2(n1.x,n2.x,w);
  }else if(method=="bayesmom"){
    float dn   = ladd(lgamma(n1.size+n2.size),n1.d+n2.d);
    float pi   = lgamma(n1.size+n2.size)-dn;
    NumericVector x   = n1.x + n2.x;
    float L    = log_dirichlet_multinom(x);
    float Lt   = ladd(pi+L,-pi+n1.Lt+n2.Lt);
    // reverse for min heap
    d = Lt-pi-L;
  }else{
    stop("Unknown method");
  } 
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
  int D = X.ncol();
  int T = sum(X);

  NumericVector w(D);
  if(method=="chi2"){
    for(int d=0; d<D; ++d){
      w(d)=T/sum(X(_,d));
    }
  }

  
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
      if(method=="bayesmom"){
        cnode.d=0;
        cnode.pi=0;
        cnode.L   = log_dirichlet_multinom(cnode.x);
        cnode.Lt  = cnode.pi+cnode.L;
        cnode.r   = 0;
      }
      for(int n=0; n<nbi.length(); ++n){
        int j = nbi[n];
        if(i!=j){
          NumericMatrix::Row x2 = X(j,_);
          node vnode;
          vnode.x = x2;
          vnode.size=1;
          if(method=="bayesmom"){
            vnode.d=0;
            vnode.pi=0;
            vnode.L   = log_dirichlet_multinom(vnode.x);
            vnode.Lt  = vnode.pi+vnode.L;
            vnode.r   = 0;
          }
          float d = dist(cnode,vnode,w,method);
          //Rcout << d << std::endl;
          cnode.neibs.insert(std::make_pair(j,d));
          if(i<j){
            priority_queue.insert(std::make_pair(d,std::make_pair(i,j)));
          }
        }
      }
      graph[i]=cnode;
    }
  }
  

  //Rcout << "--- Init end ----" << std::endl;
  // Lets Merge !
  float H =0;
  NumericMatrix merge(V-1,2);
  NumericVector height(V-1);
  
  for(int imerge=0;imerge<(V-1);imerge++){
    

    int node_id = V+imerge;
    auto best_merge = priority_queue.begin();
    //Rcout << "###############"  << std::endl;
    //print_pq(priority_queue);
    
    
    // deal with isolated regions (no more merge possible and node)
    if(best_merge==priority_queue.end()){

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
    
    
    
    Rcout << "Merge :" << g << "--" << h << std::endl;
    
    
    
    H = H + best_merge->first;
    height[imerge]=best_merge->first;
    
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
    if(method=="ward" || method=="centroid"){
        new_node.x  = (graph[g].size*graph[g].x + graph[h].size*graph[h].x)/(graph[g].size+graph[h].size);      
    }else if(method=="median"){
        new_node.x  = (graph[g].x + graph[h].x)/2;   
    }else if(method=="chi2"){
        new_node.x = graph[g].x + graph[h].x;
    }else if(method=="bayesmom"){
        new_node.d   = ladd(lgamma(graph[g].size+graph[h].size),graph[g].d+graph[h].d);
        new_node.pi  = lgamma(graph[g].size+graph[h].size)-new_node.d;
        new_node.x   = graph[g].x + graph[h].x;
        new_node.L   = log_dirichlet_multinom(new_node.x);
        new_node.Lt  = ladd(new_node.pi+new_node.L,-new_node.pi+graph[g].Lt+graph[h].Lt);
        new_node.r   = new_node.pi+new_node.L-new_node.Lt;
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
      bool found = false;
      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(v);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          priority_queue.erase(s);
          found=true;
          break;
        }
      }
      // if(!found){
      //   Rcout << "missedlink" << std::endl;
      //   Rcout << i << "-" << j << ":" << v << std::endl;
      //   Rcout << "---- pq ----" << std::endl;
      //   print_pq(priority_queue);
      // }
    

      // old link deletion in graph
      graph[j].neibs.erase(i);
      // new links in graph
      if(j!=h){
        // distance calculation
        float d = dist(new_node,graph[j],w,method);
        //Rcout << d << std::endl;
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
      
    }
    
    
    for(auto nei_h = node_h.neibs.begin();nei_h!=node_h.neibs.end();nei_h++){
      
      int i = h;
      int j = nei_h->first;
      float v = nei_h->second;
      bool found = false;

      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(v);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          priority_queue.erase(s);
          break;
        }
      }
      // if(!found){
      //   Rcout << "missedlink" << std::endl;
      //   Rcout << i << "-" << j << ":" << v << std::endl;
      //   Rcout << "---- pq ----" << std::endl;
      //   print_pq(priority_queue);
      // }
      

      // old link deletion in graph
      graph[j].neibs.erase(i);
      
      // new links in graphs
      if(j!=g){
        float d = dist(new_node,graph[j],w,method);
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
  
  // Export Centers
  NumericMatrix centers(V-1,X.ncol());
  for(int i=V;i<(2*V-1);i++){
    node cnode = graph[i];
    centers(i-V,_)=cnode.x;
  }
  CharacterVector ch = colnames(X);
  colnames(centers) = ch;
  
  return List::create(Named("merge",merge),Named("height",height),Named("data",X),Named("centers",centers));
  
}


