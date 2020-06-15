// 6/15/2020: This file is transferred from the structgen repo.

#ifndef VDWRADII_H
#define VDWRADII_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>

#include "msttypes.h"

enum atomInteraction { CLASH = 0, CONTACT = 1, INDEPENDENT = 2};

class vdwRadii {
public:
	static bool initConstants();
	static double maxSumRadii();
	static double getRadii(const string& atomName);
	static double getRadii(const Atom& a1);
	static double sumRadii(const Atom& a1, const Atom& a2);
	static bool clash(const Atom& a1, const Atom& a2, double lb = 0.7);
	static bool contact(const Atom& a1, const Atom& a2, double lb = 0.7, double ub = 1.0);
	static bool independent(const Atom& a1, const Atom& a2, double ub = 1.0);
	static atomInteraction interactionType(const Atom& a1, const Atom& a2, double lb = 0.7, double ub = 1.0);
private:
	static map<string, double> radii;
};


class resVDWRadii {
public:
	static bool initConstants();
	static double maxSumRadii();
	static double getRadii(const string& resName, const string& atomName);
	static double getRadii(Atom& a);
	static double sumRadii(Atom& a1, Atom& a2);
	static bool clash(Atom& a1, Atom& a2, double lb = 0.7);
	static bool contact(Atom& a1, Atom& a2, double lb = 0.7, double ub = 1.0);
	static bool independent(Atom& a1, Atom& a2, double ub = 1.0);
	static atomInteraction interactionType(Atom& a1, Atom& a2, double lb = 0.7, double ub = 1.0);
private:
	static map<string, map<string, double>> radii;
};

#endif
