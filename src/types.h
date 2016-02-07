
#include <iostream>
#include <cassert>
#include <string>
#include <ctime>

//STL headers
#include <vector>
#include <list>
#include <map>

//Boost headers
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

//#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/binomial_heap.hpp>

//Rcpp headers
#include <Rcpp.h>
#include <RcppEigen.h>


#define DEBUG_MODE
//#undef DEBUG_MODE


#ifdef DEBUG_MODE
#define RASSERT(condition){if(!(condition)){throw std::range_error(std::string("internal error!@")+__FILE__+":"+std::to_string((long long)__LINE__));}}
#else
#define RASSERT(condition){}
#endif


using Rcpp::Rcout;
using std::endl;
using std::ostream;
using std::string;
using std::list;
using std::vector;
using std::map;
using std::range_error;
using std::pair;

using Eigen::SparseMatrix;

//right now it supports the following linkage types
enum LinkageType{AVERAGE, SINGLE, COMPLETE};

struct EdgeObj;
struct WeightObj;

//compare functors
struct LessWeight;
struct GreaterWeight;

//NOTE: We assume we are given similarity and sets with the maximum similarity merge first
//typedef boost::heap::fibonacci_heap<WeightObj, boost::heap::compare<GreaterWeight> > MyHeap; 
typedef boost::heap::binomial_heap<WeightObj, boost::heap::compare<GreaterWeight> > MyHeap; 

typedef std::list<EdgeObj>  EdgeList;
typedef EdgeList::iterator  EdgeListItr;

typedef std::pair<int,int>  NodePair;
typedef std::map<NodePair, EdgeListItr> MyItrMap;
typedef std::map<int, int> IntIntMap;


typedef boost::numeric::ublas::matrix<int> IntegerMatrix;
//typedef boost::numeric::ublas::matrix<double> DoubleMatrix;
typedef Eigen::MatrixXd DoubleMatrix;
typedef boost::numeric::ublas::vector<int> IntegerVector;
typedef boost::numeric::ublas::vector<bool> BoolVector;

typedef std::vector<double> DoubleVector;

struct SufficientStat{
    private:
        double accumulate;
    public:
        SufficientStat():accumulate(0){}
        void operator()(double s){accumulate = s;}
        void keepMin(double s){ accumulate = accumulate > s ? s : accumulate; }
        void keepMax(double s){ accumulate = accumulate < s ? s : accumulate; }
        void keepSum(double s){ accumulate +=s;	}
        double get()const {return accumulate;}
};


struct EdgeObj{
    private:
        bool valid;
        SufficientStat stat;
    public:
        int firstNode, secondNode;
        int fNNodes, sNNodes;
        EdgeListItr shadow;

        EdgeObj(int _firstNode, int _secondNode, int _fNNodes, 
                int _sNNodes, double _eweight, EdgeListItr _shadow):
            valid(1), firstNode(_firstNode), secondNode(_secondNode), 
            fNNodes(_fNNodes), sNNodes(_sNNodes), shadow(_shadow){ 
                stat(_eweight);
            }
        EdgeObj(){}
        void set_shadow(EdgeListItr &newShadow){shadow=newShadow;}

        void do_invalid() {valid = false;}

        bool is_valid() const{return valid;}

        void update(const int newSecondNode, const int newSNNodes, EdgeListItr newShadow){
            secondNode = newSecondNode;
            shadow  = newShadow;
            sNNodes = newSNNodes;
        }

        template <LinkageType T>
            double get_weight() const{
                if(T == COMPLETE )
                    return stat.get();
                if(T == SINGLE )
                    return stat.get();
                if(T == AVERAGE )
                    return stat.get()/(fNNodes*sNNodes);
            }

        SufficientStat get_accumulator() const{ return stat; }

        template <LinkageType T>
            void accumulate_weight(SufficientStat w){ 
                if(T == AVERAGE )
                    stat.keepSum(w.get() );
                if(T == COMPLETE )
                    stat.keepMin(w.get() );
                if(T == SINGLE )
                    stat.keepMax(w.get() );
            } 

        EdgeObj get_reversed_copy(const EdgeListItr givenItr) const{ 
            EdgeObj eo = *this; 
            std::swap(eo.firstNode, eo.secondNode);
            std::swap(eo.fNNodes, eo.sNNodes);
            eo.shadow = givenItr;
            return eo;
        }
};

struct WeightObj{
    const int firstNode, secondNode;
    const EdgeListItr edgeItr1, edgeItr2;
    const double weight;
    WeightObj(const int _firstNode, const int _secondNode, 
            const EdgeListItr _edgeItr1, const EdgeListItr _edgeItr2, const double _weight):
        firstNode(_firstNode), secondNode(_secondNode), 
        edgeItr1(_edgeItr1), edgeItr2(_edgeItr2), weight(_weight){}
};

struct GreaterWeight{
    bool operator()(const WeightObj &w1, const WeightObj &w2) const
    {return w2.weight > w1.weight;}
};

struct LessWeight{
    bool operator()(const WeightObj &w1, const WeightObj &w2) const
    {return w2.weight < w1.weight;}
};

