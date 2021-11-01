#include "msttypes.h"

int main() {
  
  map<string,mstreal> test_map;
  
  // default construct a value in map and add one to it
  test_map["default"] += 1;
  
  // explicitly construct a value of 1 in map, then add 1 to it
  test_map["explicit"] = 0;
  test_map["explicit"] += 1;
  
  cout << "test_map['default']: " << test_map["default"] << endl;
  cout << "test_map['explicit']: " << test_map["explicit"] << endl;
  
  // try to check for a key before adding to it
  if (test_map.count("explicit") == 0)
  test_map["explicit"] += 1;
  
  if (test_map.count("explicit2") == 0)
  test_map["explicit2"] += 1;
  
  cout << "test_map['explicit']: " << test_map["explicit"] << endl;
  cout << "test_map['explicit2']: " << test_map["explicit2"] << endl;
  
  // seems like mstreal is default constructed to 0... is that always true?
  // https://stackoverflow.com/questions/25908281/c-double-initialization-default-values

  return 0; }
