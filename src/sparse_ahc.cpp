
#include "types.h"

// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::depends("BH")]]
// [[Rcpp::plugins(cpp11)]]


struct MapCmp {
    double eps = 0.0001;
    bool operator()(const pair<NodePair, double>& p1, const pair<NodePair, double>& p2) const{
        return (p1.second < p2.second - eps) || (p1.first != p2.first);
    }

} mCmp;

// [[Rcpp::export]]
bool dgCIsSymmetric(Eigen::SparseMatrix<double> S, double eps){

    map<NodePair, double> pool_1, pool_2;

    for (int k=0; k<S.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(S,k); it; ++it){

            if( it.row() == it.col()  || it.value() < eps )
                continue;

            if( it.row() > it.col() )
                pool_1[ NodePair(it.col(), it.row()) ] = it.value();
            else
                pool_2[ NodePair(it.row(), it.col()) ] = it.value();
        }

    std::vector<std::pair<NodePair,double> > symDifference;

    mCmp.eps = eps;
    std::set_symmetric_difference(
            pool_1.begin(), pool_1.end(),
            pool_2.begin(), pool_2.end(),
            std::back_inserter(symDifference), mCmp );

    if( symDifference.size() > 0)
        return false;

    return true;
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

        RASSERT( itrPool.find(NodePair(s,f)) != itrPool.end() );
        RASSERT( itrPool.find(NodePair(f,s)) != itrPool.end() );
        if ( f < s ) //NOTE: working with undirected graphs/no selfloops
            theHeap.push( WeightType( f, s, itrPool[NodePair(f,s)], 
                        itrPool[NodePair(s,f)], itr->get_weight<T>() ) ); 
        //setting the shadow  
        itr->set_shadow( itrPool[NodePair(s,f)] );
    }
    //print_list<T>(theList);
    return 0;
}

inline EdgeListItr go_to_begin(const EdgeList &theList, const EdgeListItr &refItr){
    auto itr = refItr;

    RASSERT(itr != theList.end());

    const int u = refItr->firstNode;
    while(itr->firstNode == u && itr != theList.begin())
        --itr;
    while( itr->firstNode != u) 
        ++itr;

    RASSERT(u == itr->firstNode );
    return itr;
}

inline EdgeListItr go_to_end(const EdgeList &theList, const EdgeListItr &refItr){
    auto itr = refItr;

    RASSERT(itr != theList.end());

    const int u = refItr->firstNode;
    while(itr != theList.end() && itr->firstNode == u )
        ++itr;
    return itr;
}

inline void fill_hierarchy_matrix(const int hCounter, const int firstNode, const int secondNode, 
        const int newNode, const double topWeight, DoubleMatrix &h){

    h(hCounter, 0) = firstNode;  
    h(hCounter, 1) = secondNode;
    h(hCounter, 2) = newNode;  			
    h(hCounter, 3) = topWeight;
}

template<LinkageType T>
int do_sparse_linkage(MyHeap &theHeap, EdgeList &theList, const int nNodes, DoubleMatrix &h){

    RASSERT(h.rows() == nNodes-1 && h.cols() == 4);

    typedef typename EdgeList::value_type  EdgeType;
    typedef typename MyHeap::value_type  WeightType;

    BoolVector legal(2*nNodes, 1);
    long newNode  = nNodes-1, hCounter = 0, garbageTime = 0;

    while( !theHeap.empty() ){/*finishing criterion*/ 

        auto fiItr = theHeap.top().edgeItr1, seItr = theHeap.top().edgeItr2;
        const int setSize = fiItr->fNNodes + fiItr->sNNodes;

        //NOTE: A paranoid check of all assumptions. Undef debug_mode macro to get ride of it.
        RASSERT( fiItr->fNNodes == seItr->sNNodes );
        RASSERT( seItr->fNNodes == fiItr->sNNodes );
        RASSERT( theHeap.top().firstNode == fiItr->firstNode ||  
                theHeap.top().secondNode == fiItr->firstNode);
        RASSERT( theHeap.top().firstNode == seItr->firstNode ||  
                theHeap.top().secondNode == seItr->firstNode);
        RASSERT( fiItr->firstNode == seItr->secondNode );
        RASSERT( fiItr->secondNode == seItr->firstNode );
        RASSERT( fiItr->is_valid() && seItr->is_valid() );


        const double topWeight = theHeap.top().weight;
        theHeap.pop(); // now discarding from the heap.
                       // invalidating it in the list too.
        fiItr->do_invalid(); 	
        seItr->do_invalid();

        ++newNode;
        fill_hierarchy_matrix(hCounter++, fiItr->firstNode, seItr->firstNode, newNode, topWeight, h);

        ///print_list<T>(theList);
        //
        for(int mode=0; mode<2; ++mode){ // now updating weights

            swap(fiItr, seItr);
            const int u=fiItr->firstNode, v=seItr->firstNode;

            map<int, EdgeList::iterator> vNodeSet;
            if( !mode ) 
                for(auto vItr = go_to_begin(theList, seItr); 
                        vItr != theList.end() && vItr->firstNode == v; ++vItr)
                    vNodeSet[vItr->secondNode] = vItr;

            for (auto uItr = go_to_begin(theList, fiItr); 
                    uItr != theList.end() && uItr->firstNode == u; ++uItr){ 

                if( !(uItr->is_valid() && legal(uItr->secondNode)) )
                    continue;
                uItr->do_invalid();

                RASSERT( uItr->secondNode != u ); //no self-loop
                RASSERT( uItr->secondNode != v ); //should not be the chosed edge

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

// Finds a permutation of tips (terminal nodes) so that the tree edges do not cross in plot.hclust.
DoubleVector order_leaves(DoubleMatrix &h, const int size){
    list<int> oList;

    map<int, list<int>::iterator> store; 

    for(int i=size-1; i>-1; --i){
        if( store.find( h(i,2) ) != store.end() ){
            auto iter     = oList.insert(store[h(i,2)], h(i,1));
            store[h(i,1)] = iter;
            iter          = oList.insert(store[h(i,2)], h(i,0));
            store[h(i,0)] = iter;
        } else {
            oList.push_front(h(i,0));
            store[h(i,0)] = oList.begin();
            oList.push_front(h(i,1));
            store[h(i,1)] = oList.begin();
        }
    }

    { //NOTE: If the input matrix represents a disconnected graph then at this step h 
      // contains multiple trees which is ok theoretically. But ``hclust'' 
      // in ``igraph'' doesn't like it and ``cutree'' ends with seg fault so here I add some
      // entries to connect these trees.
      
        vector<int> rIndices;
        for(int i=0; i<size; ++i)
            if( store.find(h(i,2)) == store.end() )
                rIndices.push_back(i);

        //NOTE: I set the hight to zero just for now!
        const double hight = 0;
        if(rIndices.size() > 1){
            int newNodes = h(size-1,2)+1, hCounter=size;
            const auto idx1 = rIndices[0]; 
            const auto idx2 = rIndices[1]; 
            fill_hierarchy_matrix(hCounter++, h(idx1,2), h(idx2,2), newNodes, hight, h);

            for(int i=2; i<rIndices.size(); ++i){
                const auto idx = rIndices[i]; 
                fill_hierarchy_matrix(hCounter++, h(idx,2), newNodes, newNodes+1, hight, h);
                ++newNodes;
            }
        }
    }

    DoubleVector ordering;
    ordering.reserve(size+1);
    for( auto iter = oList.begin(); iter != oList.end(); ++iter )
        if ( *iter < 0 )
            ordering.push_back( -1*(*iter) );

    return ordering;
}

// [[Rcpp::export]]
Rcpp::List  run_sparseAHC(
        Eigen::SparseMatrix<double> S, //similarity matrix.
        Rcpp::CharacterVector method="average"
        ){

    const int m = S.cols();
    const int n = S.rows();

    if ( n != m)
        throw range_error("input matrix must be symmetric!");

    LinkageType linkageType = AVERAGE;	
    if( Rcpp::as<string> (method) == "average" )
        linkageType = AVERAGE;
    else if( Rcpp::as<string> (method) == "single" )
        linkageType = SINGLE;
    else if( Rcpp::as<string> (method) == "complete" )
        linkageType = COMPLETE;
    else
        throw range_error("undefined method [average,single,complete]!");

    DoubleMatrix hierarchy(n-1, 4);

    const int hCounter = sparse_linkage(S, n, linkageType, hierarchy);

    //NOTE: changing hierarchy to hclust format
    for(int i=0; i<hCounter; ++i){
        for(int j=0; j<3; ++j){
            ++hierarchy(i,j); 
            hierarchy(i,j) = ( hierarchy(i,j) > n) ? hierarchy(i,j)-n : -1*hierarchy(i,j);
        }
        hierarchy(i,3) = -1*hierarchy(i,3);
    }

    auto orderedTips = order_leaves(hierarchy, hCounter);

    //Rcpp::List L =  Rcpp::List::create( Rcpp::Named("merge") = hierarchy.block(0, 0, hCounter, 2),
    //        Rcpp::Named("height") = hierarchy.block(0, 3, hCounter, 1),
    //        Rcpp::Named("method") = method,
    //        Rcpp::Named("order")  = orderedTips
    //        );
    Rcpp::List L =  Rcpp::List::create( Rcpp::Named("merge") = hierarchy.block(0, 0, n-1, 2),
            Rcpp::Named("height") = hierarchy.block(0, 3, n-1, 1),
            Rcpp::Named("method") = method,
            Rcpp::Named("order")  = orderedTips
            );
    L.attr("class")= "hclust";
    return L;
}
