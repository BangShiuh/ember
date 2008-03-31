#pragma once
#include <string>
#include "mathUtils.h"

class configOptions
{
public:
	void readOptionsFile(const std::string& filename);

	std::string inputDir;
	std::string outputDir;
	std::string restartFile;
	bool useRelativeRestartPath;

	bool overrideTu;
	bool overrideReactants;
	bool haveRestartFile;

	bool fixedBurnedVal;
	bool fixedLeftLoc;
	bool twinFlame;
	bool curvedFlame;

	int regridStepInterval;
	int outputStepInterval;
	int profileStepInterval;
	double regridTimeInterval;
	double outputTimeInterval;
	double profileTimeInterval;
	double maxTimestep;
	
	double idaRelTol;
	double idaContinuityAbsTol;
	double idaMomentumAbsTol;
	double idaEnergyAbsTol;
	double idaSpeciesAbsTol;

	std::string gasMechanismFile;
	std::string gasPhaseID;
	std::string fuel;
	std::string oxidizer;
	double equivalenceRatio;
	dvector reactants; // mole fractions
	dvector products;
	double pressure;

	double Tu, Tb;
	double strainRateInitial, strainRateFinal;
	double strainRateDt, strainRateT0;

	int nPoints;
	double xLeft, xRight;

	double vtol, dvtol, rmTol, dampConst, gridMin, gridMax;
	double uniformityTol, absvtol;
	double boundaryTol, boundaryTolRm;

	bool unburnedLeft;
	int addPointCount;

	double tStart, tEnd;

	int gridAlpha;
	int kContinuity, kMomentum, kEnergy, kSpecies;

	bool outputAuxiliaryVariables;
	bool outputTimeDerivatives;
	bool outputHeatReleaseRate;

	bool terminateForSteadyQdot;
	double terminationTolerance;
	double terminationPeriod;

	dvector strainRateList;

	int outputFileNumber; // number of output files written
	bool fileNumberOverride; // true if outputFileNumbe was given in the input file
};
