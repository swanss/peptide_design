#ifndef ranklist_h
#define ranklist_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <list>

/**
 Utility class that maintains a list of k max (or min) elements ranked by some
 key.
 */
template <class T>
class RankList {
public:
    /**
     Initialize a rank list.
     
     @param size the number of elements to hold
     @param maximize if true, store the items with max scores; otherwise, store
            the mins
     */
    RankList(int size = 1, bool maximize = true): _size(size), _max(maximize) {};

    /// @return the current best items along with their scores
    vector<pair<T, double>> values() {
        vector<pair<T, double>> ret;
        ret.insert(ret.end(), _values.begin(), _values.end());
        return ret;
    }
    
    /// @return the current best items without their scores
    vector<T> items() {
        vector<T> result;
        //transform(_values.begin(), _values.end(), result.begin(), [](pair<T, double> value) { return value.first; });
        for (pair<T, double> val: _values) {
            result.push_back(val.first);
        }
        return result;
    }
    
    /**
     Insert the given item into the rank list if its score is sufficiently good.
     
     @param item the item to potentially insert
     @param score the score by which to compare the item to those in the rank list
     @return true if the item was inserted, and false if it did not make the
             rankings
     */
    bool insert(T item, double score) {
        auto insertionIt = _values.end();
        for (auto it = _values.begin(); it != _values.end(); ++it) {
            if ((_max && score > (*it).second) || (!_max && score < (*it).second)) {
                insertionIt = it;
                break;
            }
        }
        
        if (insertionIt == _values.end() && _values.size() >= _size)
            return false;
        
        _values.insert(insertionIt, pair<T, double>(item, score));
        if (_values.size() > _size)
            _values.pop_back();
        
        return true;
    }
    
    /**
     Remove and return the best item in the rank list. Items that do not make
     the ranking are not saved, so the size of the rank list will decrease by
     one.
     
     @return the best element in the rank list
     */
    T pop() {
        T val = _values.front().first;
        _values.pop_front();
        return val;
    }
    
    /**
     Subscript operator
     */
    T& operator[] (const int index) {
        auto it = _values.begin();
        for (int i = 0; i < index; i++)
            ++it;
        return (*it).first;
    }
    
    /**
     @return the current size of the rank list
     */
    int size() {
        return _values.size();
    }
    
private:
    int _size = 10;
    bool _max = true;
    
    list<pair<T, double>> _values;
};


#endif /* ranklist_h */
