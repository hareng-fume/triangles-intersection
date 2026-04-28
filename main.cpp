#include "triangle_intersect.h"

#include <iostream>
#include <vector>

using namespace triangle_3d;

struct TestCase
{
  int index;
  bool expected;
  bool actual;

  double t1[9];
  double t2[9];

  void Run() { actual = intersect(t1, t2); }
  bool Success() const {return expected == actual;}
};


int main()
{
  std::vector<TestCase> test_case = {
    {0, true, false, {0,0,0,2,0,0,0,2,0}, {0.5,0.5,-1,0.5,0.5,1,1.5,0.5,0}},      // 1.  Plane crossing center
    {1, true, false, {0,0,0,2,0,0,0,2,0}, {1,-1,-1,1, 1, 1,3, 0, 0}},             // 2.  Edge crossing plane
    {2, true, false, {0,0,0,2,0,0,0,2,0}, {1,1,0,1,1,1,2,2,2}},                   // 3.  Vertex touches triangle plane
    {3, false, true, {0,0,0,2,0,0,0,2,0}, {0,0,5,2,0,5,0,2,5}},                   // 4.  Clearly separated (parallel planes) 
    {4, false, true, {0,0,0,2,0,0,0,2,0}, {10,10,0,12,10,0,10,12,0}},             // 5.  Parallel but offset (no overlap)
    {5, true, false, {0,0,0,3,0,0,0,3,0}, {1,1,0,4,1,0,1,4,0}},                   // 6.  Coplanar overlapping triangles
    {6, false, true, {0,0,0,1,0,0,0,1,0}, {2,2,0,3,2,0,2,3,0}},                   // 7.  Coplanar but disjoint
    {7, true, false, {0,0,0,2,0,0,0,2,0}, {2,0,0,2,2,0,4,0,0}},                   // 8.  Edge-touching coplanar
    {8, true, false, {0,0,0,2,0,0,0,2,0}, {2,0,0,3,1,0,1,3,0}},                   // 9.  Vertex-touching coplanar
    {9, true, false, {0,0,0,1,0,0,0,1,0}, {0.5,0.5, 1e-10,0.5,0.5,-1e-10,2,0,0}}, // 10. Nearly degenerate (EPS stress test)
    {10,false, true, {0,0,0,2,0,0,0,2,0}, {0.5,0.5,0.01,0.5,0.5,0.02,2,0,0}},     // 11. Thin miss (very important bug finder) [failed]
    {11,true, false, {1,2,3,4,5,6,7,2,1}, {3,3,2,2,1,5,6,4,0}}                    // 12. Random skewed intersection
  };

  for(std::size_t i =0; i < test_case.size(); ++i)
  {
    test_case[i].Run();
    if(!test_case[i].Success())
	  {
      std::cout << std::boolalpha
        << "test_" << test_case[i].index
        << " failed (expected: " << test_case[i].expected
        <<          ", actual: " << test_case[i].actual
        << ")\n";
    }
  }

  return 0;
}
