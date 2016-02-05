
#include "types.h"

// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::depends("BH")]]
// [[Rcpp::plugins(cpp11)]]

template <LinkageType T>
inline void print_list(const EdgeList &theList)
{
    Rcout<<"+>-------------------\n";
    for(auto itr= theList.begin(); itr != theList.end(); ++itr)
    {
        Rcout << itr->firstNode << "\t";
        Rcout << itr->secondNode << "\t";
        Rcout << itr->get_weight<T>() << "\t";
        Rcout << itr->shadow->firstNode << "\t";
        Rcout << itr->shadow->secondNode << "\t";
        Rcout << itr->is_valid() << "\t";
        Rcout<<endl;
    }
    Rcout<<"-<-------------------\n";
} 

template <LinkageType T>
bool fill_heap_and_list(SparseMatrix<double> S, const int nNodes, MyHeap &theHeap, EdgeList &theList){

    typedef typename EdgeList::value_type EdgeType;
    typedef typename MyHeap::value_type WeightType;


    MyItrMap itrPool;

    for (int k=0; k<S.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(S,k); it; ++it){
            if( it.row() == it.col() ) //NOTE: no selfloop !?
                continue;
            //endpoint1,  endpoint2, nNodes in endpoint1, nNodes in endpoint2, weight, shadow(will be set later)
            theList.push_back( EdgeType(it.col(), it.row(), 1, 1, it.value(), theList.begin()) );
            itrPool[ NodePair(it.col(), it.row()) ]  = --theList.end();
        }

    auto itr = theList.begin();
    for (; itr != theList.end(); ++itr){
        const int f = itr->firstNode, s = itr->secondNode;

#ifdef DEBUG_MODE
        RASSERT( itrPool.find(NodePair(s,f)) != itrPool.end() );
        RASSERT( itrPool.find(NodePair(f,s)) != itrPool.end() );
#endif
        if ( f < s ) //NOTE: working with undirected graphs/no selfloops
            theHeap.push( WeightType( f, s, itrPool[NodePair(f,s)], 
                        itrPool[NodePair(s,f)], itr->get_weight<T>() ) ); 
        //setting the shadow  
        itr->set_shadow( itrPool[NodePair(s,f)] );
    }

#ifdef DEBUG_MODE
    //print_list<T>(theList);
#endif
    return 0;
}

inline EdgeListItr go_to_begin(const EdgeList &theList, const EdgeListItr &refItr){
    auto itr = refItr;
#ifdef DEBUG_MODE
    RASSERT(itr != theList.end());
#endif
    const int u = refItr->firstNode;
    while(itr->firstNode == u && itr != theList.begin())
        --itr;
    while( itr->firstNode != u) 
        ++itr;
#ifdef DEBUG_MODE
    RASSERT(u == itr->firstNode );
#endif
    return itr;
}

inline EdgeListItr go_to_end(const EdgeList &theList, const EdgeListItr &refItr){
    auto itr = refItr;
#ifdef DEBUG_MODE
    RASSERT(itr != theList.end());
#endif
    const int u = refItr->firstNode;
    while(itr != theList.end() && itr->firstNode == u )
        ++itr;
    return itr;
}

inline void fill_hierarchy_matrix(const int hCounter, const int firstNode, const int secondNode, 
        const int newNode, const double topWeight, DoubleMatrix &h){
#ifdef DEBUG_MODE
    //cout <<"("<<firstNode<<","<<secondNode<<","<<topWeight<<") --> "<<newNode<<endl;
#endif
    h(hCounter, 0) = firstNode;  
    h(hCounter, 1) = secondNode;
    h(hCounter, 2) = newNode;  			
    h(hCounter, 3) = topWeight;
}

template<LinkageType T>
int do_sparse_linkage(MyHeap &theHeap, EdgeList &theList, const int nNodes, DoubleMatrix &h){

#ifdef DEBUG_MODE
    RASSERT(h.rows() == nNodes-1 && h.cols() == 4);
#endif

    typedef typename EdgeList::value_type  EdgeType;
    typedef typename MyHeap::value_type  WeightType;

    BoolVector legal(2*nNodes, 1);
    long newNode  = nNodes-1, hCounter = 0, garbageTime = 0;

    while( !theHeap.empty() ){/*finishing criterion*/ 

        auto fiItr = theHeap.top().edgeItr1, seItr = theHeap.top().edgeItr2;
        const int setSize = fiItr->fNNodes + fiItr->sNNodes;

#ifdef DEBUG_MODE
        RASSERT( fiItr->fNNodes == seItr->sNNodes );
        RASSERT( seItr->fNNodes == fiItr->sNNodes );
        RASSERT( theHeap.top().firstNode == fiItr->firstNode ||  
                theHeap.top().secondNode == fiItr->firstNode);
        RASSERT( theHeap.top().firstNode == seItr->firstNode ||  
                theHeap.top().secondNode == seItr->firstNode);
        RASSERT( fiItr->firstNode == seItr->secondNode );
        RASSERT( fiItr->secondNode == seItr->firstNode );
        RASSERT( fiItr->is_valid() && seItr->is_valid() );
#endif
        const double topWeight = theHeap.top().weight;
        theHeap.pop(); // now discarding from the heap.
                       // invadiating it in the list too.
        fiItr->do_invalid(); 	
        seItr->do_invalid();

        ++newNode;
        fill_hierarchy_matrix(hCounter++, fiItr->firstNode, seItr->firstNode, newNode, topWeight, h);
#ifdef DEBUG_MODE
        ///print_list<T>(theList);
#endif
        for(int mode=0; mode<2; ++mode){ // now updating weights

            swap(fiItr, seItr);
            const int u=fiItr->firstNode, v=seItr->firstNode;

            std::map<int, EdgeList::iterator> vNodeSet;
            if( !mode ) 
                for(auto vItr = go_to_begin(theList, seItr); 
                        vItr != theList.end() && vItr->firstNode == v; ++vItr)
                    vNodeSet[vItr->secondNode] = vItr;

            for (auto uItr = go_to_begin(theList, fiItr); 
                    uItr != theList.end() && uItr->firstNode == u; ++uItr){ 

                if( !(uItr->is_valid() && legal(uItr->secondNode)) )
                    continue;
                uItr->do_invalid();
#ifdef DEBUG_MODE
                RASSERT( uItr->secondNode != u ); //no self-loop
                RASSERT( uItr->secondNode != v ); //should not be the chosed edge
#endif
                const int neiNode  = uItr->secondNode; 
                const auto connTie = uItr->shadow;

                if ( !connTie->is_valid() ) 
                    continue;

                if( !mode  && vNodeSet.find( neiNode ) != vNodeSet.end() ){
                    auto vItr = vNodeSet[ neiNode]; 
                    if( vItr->is_valid() )
                        connTie->accumulate_weight<T>(vItr->get_accumulator());
                    vItr->do_invalid();
                }

                //NOTE: adding a dummy obj	
                theList.push_back( EdgeType() );	
                connTie->update(newNode, setSize, --theList.end());
                theList.back() = connTie->get_reversed_copy(connTie);
                theHeap.push( WeightType( connTie->firstNode, connTie->secondNode, 
                            --theList.end(), connTie, connTie->get_weight<T>()) ); 
            }
        }

        //NOTE: cleaning the list to save memory!
        theList.erase( go_to_begin(theList, fiItr), go_to_end(theList, fiItr) ); 
        theList.erase( go_to_begin(theList, seItr), go_to_end(theList, seItr) ); 
        legal(fiItr->firstNode )  = legal(seItr->firstNode)  = 0;

        //NOTE: throwing away invalid weights at the top of the heap
        while( !theHeap.empty()  )
            if( legal(theHeap.top().firstNode) && legal(theHeap.top().secondNode) 
                    && theHeap.top().edgeItr1->is_valid() && theHeap.top().edgeItr2->is_valid() )
                break;
            else
                theHeap.pop();

        //NOTE: cleaning the heap to save some memory
        // can we do this part in parallel to the previous part?
        // we can do that whenever it is necessary too but I don't think it takes too long.
        if ( ++garbageTime % 2000 == 0){
            MyHeap tmpHeap;
            while( !theHeap.empty()  ){
                if( legal(theHeap.top().firstNode) && legal(theHeap.top().secondNode) 
                        && theHeap.top().edgeItr1->is_valid() && theHeap.top().edgeItr2->is_valid() )
                    tmpHeap.push( theHeap.top() );
                theHeap.pop();
            }

            while( !tmpHeap.empty() ){
                theHeap.push( tmpHeap.top() );
                tmpHeap.pop();	
            }
            garbageTime = 0;
        } //end of garbageTime
    }//end of main while
    return hCounter;
}

int sparse_linkage(SparseMatrix<double> S, const int nNodes, const LinkageType t, DoubleMatrix &h){

#ifdef DEBUG_MODE
    Rcout<<"warning: running in debug mode!" <<endl;
#endif
    MyHeap theHeap;
    EdgeList theList;

    switch(t){
        case AVERAGE:
            fill_heap_and_list<AVERAGE>(S, nNodes, theHeap, theList);
            return do_sparse_linkage<AVERAGE>(theHeap, theList, nNodes, h);

        case COMPLETE:
            fill_heap_and_list<COMPLETE>(S, nNodes, theHeap, theList);
            return do_sparse_linkage<COMPLETE>(theHeap, theList, nNodes, h);

        case SINGLE:
            fill_heap_and_list<SINGLE>(S, nNodes, theHeap, theList);
            return do_sparse_linkage<SINGLE>(theHeap, theList, nNodes, h);

        default:
            RASSERT(false);
    }
}

DoubleVector order_leaves(DoubleMatrix &h, const int size){
    std::list<int> oList;
    std::map<int, std::list<int>::iterator> store; 

    for(int i=size-1; i>-1; --i){

        if( store.find( h(i,2) ) != store.end() ){
            auto iter    = oList.insert(store[h(i,2)], h(i,1));
            store[h(i,1)]= iter;
            iter    = oList.insert(store[h(i,2)], h(i,0));
            store[h(i,0)]= iter;
        } else {
            oList.push_front(h(i,0));
            store[h(i,0)] = oList.begin();
            oList.push_front(h(i,1));
            store[h(i,1)] = oList.begin();
        }
    }

    std::vector<double> ordering;
    ordering.reserve(size);
    for( auto iter = oList.begin(); iter != oList.end(); ++iter )
        if ( *iter < 0 )
            ordering.push_back( -1*(*iter) );
    return ordering;
}

// [[Rcpp::export]]
Rcpp::List  run_sparseAHC(
        Eigen::SparseMatrix<double> S, //similarity matrix.
        Rcpp::CharacterVector method="average",
        bool noOrder=false //for huge input you may trun it on
        ){

    const int m = S.cols();
    const int n = S.rows();

    if ( n != m)
        throw std::range_error("input matrix must be symmetric!");

    LinkageType linkageType = AVERAGE;	
    if( Rcpp::as<std::string> (method) == "average" )
        linkageType = AVERAGE;
    else if( Rcpp::as<std::string> (method) == "single" )
        linkageType = SINGLE;
    else if( Rcpp::as<std::string> (method) == "complete" )
        linkageType = COMPLETE;
    else
        throw std::range_error("undefined method [average,single,complete]!");

    DoubleMatrix hierarchy(n-1, 4);

    const int hCounter = sparse_linkage(S, n, linkageType, hierarchy);


    //clock_t begin = std::clock();
    //clock_t end   = std::clock();
    //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //Rcout<<elapsed_secs<<endl;

    //NOTE: changing hierarchy to hclust format
    for(int i=0; i<hCounter; ++i){
        for(int j=0; j<3; ++j){
            ++hierarchy(i,j); 
            hierarchy(i,j) = ( hierarchy(i,j) > n) ? hierarchy(i,j)-n : -1*hierarchy(i,j);
        }
        hierarchy(i,3) = -1*hierarchy(i,3);
    }

    Rcpp::List L =  Rcpp::List::create( Rcpp::Named("merge") = hierarchy.block(0, 0, hCounter, 2),
            Rcpp::Named("height") = hierarchy.block(0, 3, hCounter, 1),
            Rcpp::Named("method") = method,
            Rcpp::Named("order")  = noOrder?DoubleVector():order_leaves(hierarchy, hCounter) 
            );
    L.attr("class")= "hclust";
    return L;
}
