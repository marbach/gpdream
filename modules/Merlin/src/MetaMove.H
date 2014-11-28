/*The goal of this class is to describe a graph edit move we make to generate a new graph.
 * We will make a move on a set of variables, typically a pair of variables, but may be more
 * if we are doing exact search. We need to store the type of move and the score
 * because of this new move. The move could be an edge addition, edge deletion, edge reversal,
 * Markov blanket addition, which is equal to addition of a set of edges */
#ifndef _MOVE_
#define _MOVE_

#include "CommonTypes.H"
class Potential;
class MetaMove
{
	public:
	MetaMove();
	~MetaMove();
	int setScoreImprovement(double);
	int setSrcMBScore(double);
	int setTargetMBScore(double);
	int setScore(double);
	int setSrcVertex(int);
	int setTargetVertex(int);
	int setConditionSet(INTINTMAP&);
	int setConditionSetInd(int);
	int setDestPot(Potential*);

	int getSrcVertex();
	int getTargetVertex();
	INTINTMAP& getConditionSet();
	int getConditionSetInd();
	double getScore();
	double getSrcMBScore();
	double getTargetMBScore();
	double getScoreImprovement();
	Potential* getDestPot();
	INTDBLMAP pll;
	private:
	double scoreDelta;
	double mbscore;
	double targetMBScore;
	int src;
	int target;
	INTINTMAP conditionSet;
	int conditionSetInd;
	Potential* destPot;
};

#endif
