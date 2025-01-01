/*
icpc *.cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -O
*/

#include "../crysPara.hpp"

int main(int argc, char** argv) {

    //start MPI    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	
    //compute resolution of each unique structure factor, output uniqStruFactReso
    computeUniqStruFactReso();
    
    /***************************************************************************

                compute true phases for all allowed origin translation

    ***************************************************************************/
    //read in Fmodel got from Phenix (and CCP4) with bulk solvent correction, output struFactAmplCal and struFactPhase
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_fmodel.hkl,.HKL,.sca";
    inputFmodel(fileNameStream.str());

    //generate true phases for all allowed origin trans, input uniqStruFactAmplCal and uniqStruFactPhase, output truePhaseOrigTran
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << PDB_CODE << "_true_phase_orig_tran_"; //.txt
    generateTruePhaseForAllOrigTran("not output true phases", fileNameStream.str());


    /***************************************************************************

                compute true prot mask for all allowed origin translation

    ***************************************************************************/

    //read in Fmodel got from Phenix (and CCP4) with bulk solvent correction, output struFactAmplCal and struFactPhase
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_fatom.hkl,.HKL,.sca";
    inputFmodel(fileNameStream.str());

    //generate true prot mask for all origin trans, input uniqStruFactAmplCal and uniqStruFactPhase, output trueProtMaskOrigTran
    generateTrueProtMaskForAllOrigTran();
    

    /***************************************************************************

                input observed structure factor amplitude

    ***************************************************************************/
    //read in observerd diffraction data in experiment, reserved in uniqStruFactStat, uniqStruFactAmplObs, uniqStruFactAmplObsSigm
    //data format (3%5d, %2d, %9.1f, %7.1f, %8.3f) for h, k, l, obsStatus, obssf, obssfSigma, resolution
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_uniq_sf.txt";
    inputUniqStruFactAmplObs(fileNameStream.str());


    /***************************************************************************

                input observed reference histogram

    ***************************************************************************/
    
    //read in the boundary values of reference protein hist, reserved in stdBoundLargeBin
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE_REFE_HIST << "_hist.txt";
    inputRefeHist(fileNameStream.str());

    /***************************************************************************

                input SigmWeigObs

    ***************************************************************************/
    /*
    //read in the sigmweigobs
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_iter_sigm.txt";
    inputSigmWeigObsAndSigmWeigAvg(fileNameStream.str());
    */
    /***************************************************************************

                initialize random density in asu

    ***************************************************************************/
    /*
    //read in density in asu from deposited pdb file, each atom is represented by density 1.0
    unsigned long seed = mixClockTimeGetid();
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_rand_seed_init_dens_asu.txt";
    ofstream outputSeed(fileNameStream.str().c_str());
    outputSeed << seed << endl;    outputSeed.close();
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << ".pdb";
    srand (seed); //initialize random seed
    inputDensAsuFromPdb(fileNameStream.str(), 0.005); //input random 2% atom positions
    */
    /*
    //generate random density in asu
    unsigned long seed = mixClockTimeGetid();
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_rand_seed_init_dens_asu.txt";
    ofstream outputSeed(fileNameStream.str().c_str());
    outputSeed << seed << endl;    outputSeed.close();
    generateRandDensAsu(seed);
    */
    //read in the boundary values of reference protein hist, reserved in stdBoundLargeBin
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_rand_seed_init_dens_asu.txt";
    unsigned long seed = inputRandomnumber(fileNameStream.str());
    generateRandDensAsu(seed);
    
    //apply symmetric operation to get density in unit cell from densAsu, input densAsu, output densUnitCell
    applySymmOperDensAsu();

    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply backward fft, input densUnitCell, output allStruFactAmpl and allStruFactPhase
    applyFFTDensToStruFact();

    //catch uniqStruFactAmplCal and uniqStruFactPhase from allStruFactAmpl and allStruFactPhase
    catchUniqStruFactAmplCalAndPhase();
    
    /***************************************************************************

                start iteration on density and phase operation

    ***************************************************************************/
    
    for(int iter = 0; iter <= numIter - 1; iter++) {
        
        //vary solvent content, input solvContIni, solvContFina, output new solvCont
        solvCont = varySolvContLinearly(iter);

        //vary sigmWeigAvg
        sigmWeigAvg = varySigmWeigAvgLinearly(iter);
        
        //vary sigmaWeigObs
        //sigmWeigObs = varySigmWeigObsNonLinearly(iter);

        //vary Fobs, input uniqStruFactAmplOrigObs, output uniqStruFactAmplObs
        //varyWeightUniqStruFactAmplObs(sigmWeigObs);
        
        //fill up missing reflections, input uniqStruFactAmplCal and uniqStruFactAmplObs, output missStruFactAmpl
        fillMissRefl(iter);

        //compute allStruFactAmpl from uniqStruFactAmplObs(and missStruFactAmpl)
        computeAllStruFactAmplFromUniqStruFactAmplObs();
        
        //compute allStruFactPhase from uniqStruFactPhase
        computeAllStruFactPhaseFromUniqStruFactPhase();

        //apply forward fft, input allStruFactAmpl and allStruFactPhase, output densUnitCell
        applyFFTStruFactToDens();

        //get densAsuTemp from densUnitCell, input densUnitCell, output densAsuTemp
        catchDensAsuTempFromDensUnitCell();
        
        //compute R values, input densAsuTemp with protein mask constraint, output Rwork
        computeRworkAndRfree(iter);
        
        //use updated deltaRwork to check the num of elites on all ranks
        updateDeltaRworkAndNumOfElite(iter);
        
        //early stop when all densAsu converge
        if(numElite == commSize && iter < numIter-numIterLimitDens-numIterSolvFlat)
            iter = numIter-numIterLimitDens-numIterSolvFlat; 
        
        //apply genetic algorithm to improve densAsuTemp, input densAsuTemp, output new densAsuTemp
        if(iter % 100 == 0 && iter > 0 && iter <= numIter - numIterSolvFlat){
            geneticAlgorithm(iter); 
        }
        
        //compute weighted average density in unit cell, input densUnitCell, output weigAvgDensUnitCell
        computeWeigAvgDensUnitCell(sigmWeigAvg);

        //get weigAvgDensAsu from weigAvgDensUnitCell, input weigAvgDensUnitCell, output weigAvgDensAsu
        catchWeigAvgDensAsuFromWeigAvgDensUnitCell();

        //compute a densCutoff value on weigAvgDensAsu by bins, input weigAvgDensAsu, return densCutoff
        densCutoff = computeWeigAvgDensCutoffAsu(1.0-solvCont);

        //make a protein mask in asu, input weigAvgDensAsu, output protMaskAsu
        makeProtMaskAsu(densCutoff);
        
        //apply HIO in asu to update densAsu, input previous densAsu and current densAsuTemp, output new densAsu
        applyHIO(iter); //HIODensLimit
        
        //apply histogram matching on densAsu, input densAsu, output new densAsu
        applyHistMatch();
        
        //apply solvent flattening on densAsu during the last iterations, input densAsu, output new densAsu
        if(iter > numIter - numIterSolvFlat) applySolvFlat();

        //apply symmetric operation to get densUniCell from updated densAsu, input densAsu, output densUnitCell
        applySymmOperDensAsu();

        //apply backward fft, input densUnitCell, output allStruFactAmpl and allStruFactPhase
        applyFFTDensToStruFact();

        //catch uniqStruFactAmplCal and uniqStruFactPhase from allStruFactAmpl and allStruFactPhase
        catchUniqStruFactAmplCalAndPhase();

        //compute mean phase error for all allowed origin translations in resolution sphere
        float meanPhaseErrorResoSphere[numResoShell][numOrigTran] = {};
        computeMeanPhaseErrorResoSphere(meanPhaseErrorResoSphere);
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell){
            for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
                meanPhaseErrorOrigTranResoSphere[iter][ishell][iorig]= meanPhaseErrorResoSphere[ishell][iorig];
        }

        //compute mean phase error for all allowed origin translations in resolution shell
        float meanPhaseErrorResoShell[numResoShell][numOrigTran] = {};
        computeMeanPhaseErrorResoShell(meanPhaseErrorResoShell);
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell){
            for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
                meanPhaseErrorOrigTranResoShell[iter][ishell][iorig]= meanPhaseErrorResoShell[ishell][iorig];
        }
        
        //compute correlation coefficient for all allowed origin translations
        float corrCoef[numResoShell][numOrigTran] = {};
        computeCorrCoefInSphere(corrCoef);
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell){
            for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
                corrCoefOrigTran[iter][ishell][iorig]= corrCoef[ishell][iorig];
        }
        
        //compute prot mask error for all allowed origin translations
        float protMaskMatch[numOrigTran] = {};
        computeProtMaskMatch(protMaskMatch);
        for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
            protMaskMatchOrigTran[iter][iorig]= protMaskMatch[iorig];

        //check origin choice
        origChoice[iter] = checkOrigChoiceByProtMaskMatch(protMaskMatch);
        
        //check convergence
        checkConvByMeanPhaseError(meanPhaseErrorResoSphere); //actually by meanPhaseErrorResoSphere[0][iorig]
        
        //compute the mean phase error of converged densAsuTemp
        computeMeanPhaseErrorOfConvDensAsuTemp(iter); //time consuming
        
        //output the progress
        if(iter % 10 == 0) cout << commRank << " / " << iter << " / " << numIter << endl;
        
    }

    /***************************************************************************

                output calculated results

    ***************************************************************************/
    //output R vlaue
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_R_work.txt";
    outputRvalue(fileNameStream.str(), Rwork);
    
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_R_free.txt";
    outputRvalue(fileNameStream.str(), Rfree);
     
    //--------------------------------------------------------------------------
    //output protein mask match
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_prot_mask_match.txt";
    outputProtMaskMatch(fileNameStream.str(), 0);
    
    //--------------------------------------------------------------------------
    //output mean phase error in resolution sphere
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_mean_phase_error_reso_sphere.txt";
    outputMeanPhaseErrorResoSphere(fileNameStream.str(), 0);
    
    //--------------------------------------------------------------------------
    //output mean phase error in resolution shell
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_mean_phase_error_reso_shell.txt";
    outputMeanPhaseErrorResoShell(fileNameStream.str(), 0);
    
    //--------------------------------------------------------------------------
    //output correlation coefficient
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_corr_coef.txt";
    outputCorrCoef(fileNameStream.str(), 0);
     
    //--------------------------------------------------------------------------
    //output mean phase error in resolution sphere of converged densAsuTemp_avg
    ofstream outputPhaseErrorResoSphereOfConvDensAsuTemp;
    if(commRank == 0) {
        fileNameStream.str( string() );
        fileNameStream.clear();
        fileNameStream << commRank << "_mean_phase_error_reso_sphere_of_conv_dens.txt";
        outputMeanPhaseErrorOfConvDensAsuTemp(fileNameStream.str(), 0);
    }
    /*
    //--------------------------------------------------------------------------
    //output the final protMask in unit cell to a pdb file
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_prot_mask_fina.pdb";
    densCutoff = computeWeigAvgDensCutoffAsu(1.0-solvCont);
    outputWeigAvgDensUnitCell(fileNameStream.str(), densCutoff);

    //--------------------------------------------------------------------------
    //output the final density in unit cell to a pdb file
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_densUnitCell_fina.pdb";
    densCutoff = computeDensCutoffUnitCell((1.0-solvCont)*0.05);
    outputDensUnitCell(fileNameStream.str(), densCutoff);
    */
    //--------------------------------------------------------------------------
    //output uniqStruFactPhase and missStruFactAmpl
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << commRank << "_" << PDB_CODE << "_phase.cif";
    outputCifFileForCoot(fileNameStream.str());
    
    //--------------------------------------------------------------------------
    //finalize MPI
    MPI_Finalize();
    
    return 0;
}
