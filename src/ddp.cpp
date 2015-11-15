#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP DDP(SEXP yvec, SEXP Xmat, SEXP Zmat,  SEXP personsvec, SEXP Dosemat, 
        SEXP numTreat, SEXP typeTreat, SEXP OmegaList, SEXP omegapluslist,
        SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP shapealph, SEXP ratebeta, SEXP Minit)
{
BEGIN_RCPP
    NumericVector yr(yvec);
    NumericMatrix Xr(Xmat);
    NumericMatrix Zr(Zmat);
    NumericMatrix Dr(Dosemat); /* np x ntr */
    List Omega(OmegaList); /* List of Omega matrices for CAR priors */
    List omegaplus(omegapluslist);
    IntegerVector numT(numTreat); /* vector of number of each treatment type, including intercept */
    // vector of covariance form choices for each treatment type under a multivariate Gaussian.
    // 1 = CAR, 2 = MVN, 3 = IND
    IntegerVector typeT(typeTreat); /* same length as numT */
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    int M     = as<int>(Minit);
    int nkeep = (niter -  nburn)/nthin;
    int nOmega = Omega.size();
    double ac = as<double>(shapealph);
    double bc = as<double>(ratebeta);

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), nty = typeT.size();
    int ntr = sum(numT) + 1; /* the added 1 for an intercept */
    /* nty = number of treatment types; ntr = total number of dosages, including intercept */
    
    // Check for entry / dimension on number of treatments
    if( typeT.size() != numT.size() )
        throw std::range_error("Number of treatments must equal number of treatment types");

    // Compute np, number of unique persons
    /* algorithm assumes cases clustered by unique client */ 
    IntegerVector dpr = diff(pr);
    icolvec diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    
    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    mat dosemat(Dr.begin(), np, ntr, false);
    colvec y(yr.begin(), nc, false);
    icolvec persons(pr.begin(), nc, false);
    icolvec numt(numT.begin(), nty, false); /* used nty since numT excludes intercept */
    icolvec typet(typeT.begin(), nty, false);
    field<mat> doseperson(np,1); /* nr x (nr*ntr) dosage matrix, by person */
    buildD(doseperson, dosemat, nr); /* build doseperson, nr x (nr*ntr) matrices for each of np persons */
    // set up data objects indexed by clustering/subject unit
    field<mat> zsplit(np,1); 
    field<mat> xsplit(np,1); 
    field<colvec> ysplit(np,1); field<colvec> ytilsplit(np,1);
    field<mat> dosequad(np,1), zbydose(np,1);
    CovUnitSplitDDP(zmat, xmat, doseperson, y, zsplit, xsplit, ysplit, dosequad, zbydose, persons);
    mat xtx; prodMatOne(xmat,xtx); 
    // capture number of treatment types allocated to each covariance matrix type.
    // use to dimension return objects
    uvec poscar = find(typet == 1); int numcar = poscar.n_elem;
    uvec posmvn = find(typet == 2); int nummvn = posmvn.n_elem;
    uvec posind = find(typet == 3); int numind = posind.n_elem;

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    //arma_rng::set_seed_random(); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    colvec Alpha(nkeep);
    colvec Conc(nkeep);
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    mat Beta(nkeep,nf);
    // random effect sampling and return structures 
    int i, dim;
    field<mat> DoseEffects((nty+1),1); /* List element is typet(m).  nkeep x np*(nr*dim_m) matrix of sampled random effects */
    for(i = 0; i < nty; i++)
    {
        dim = numt(i);
        DoseEffects(i,0).set_size(nkeep,np*nr*dim);
    }
    DoseEffects(nty,0).set_size(nkeep,nr*np); /* intercept - eqivalent to B in DP model */
    field<mat> dmats(np,1); /* nr x ntr dose effects matrix for each of np persons */
    field<mat> dosemats(np,nty); /* nr x dim_m effects matrix for trt type m for each of np persons */
    field<colvec> dosevecs(np,nty); /* nr*dim_m x 1 effects vector for trt type m for each of np persons */
    rowvec doserow; /* 1 x np*(nr*dim_m)) to capture dose parms for type(m) on a given MCMC sample */
    /* intercept random effects - equivalent to b_{1}, .., b_{q} under DP model */ 
    rowvec bint(np*nr); bint.zeros();
    /* thetamat is the np x nr person effect parameter matrix obtained by multiplying the nr x (nr*ntr)
       design matrix, doseperson(i,0) by a (nr*ntr) x 1 member of dmat (np x (nr*ntr)). */
    mat Theta(nkeep, nr*np);
    rowvec thetarow(np*nr); thetarow.zeros();
    
    // Precision parameters sampling and return objects
    /* sampling objects */
    field<mat> pmatsmvn; /* mvn */
    field<rowvec> taus; /* ind */
    field<mat> omegamats;  /* CAR */
    field<rowvec> omegapluses; /* CAR */
    rowvec alphacar; /* CAR */
    rowvec taucar; /* CAR */
    colvec logf; /* CAR - alphacar */
    field<mat> qmats; /* CAR - alphacar */
    mat pmvn; /* mvn */
    /* return objects */
    field<mat> Pmvn; /* mvn */
    field<mat> Tauind; /* ind */
    mat Alphacar; /* CAR */
    mat Taucar; /* CAR */
    field<mat> CARprecision(2,1); /* CAR */
    // set initial values for sampling objects and dimensions for return objects 
    if(nummvn > 0) /* return precision matrix objects for mvn structures */
    {
       //Define field object of mvn covariance matrices to sample */
       pmatsmvn.set_size(nummvn,1); /* sampling object */
       Pmvn.set_size(nummvn,1);   /* return object */ 
       // fill in each field member with a vector of tau values for that type */   
        for(i = 0; i < nummvn; i++)
        {
            int dosei = numt(posmvn(i)); /* number of treatments for type i */
	    Pmvn(i,0).set_size(nkeep,(dosei * dosei));
            // Pmvn(i,0).set_size(nkeep,pow(dosei,2));
            pmvn = 0.1*eye<mat>(dosei,dosei);
            pmatsmvn(i,0) = pmvn;    /* initialize sampling object */
        } 
    }
    
    if(numcar > 0) /* return CAR parameters for each CAR structure */
    {
        /* set sizes for adjacency matrix (arma) structures */
        if( numcar != nOmega )
                throw std::range_error("CAR entries not equal to number of Omega inputs");
        omegamats.set_size(numcar,1); 
        omegapluses.set_size(numcar,1); 
        for(i = 0; i < nOmega; i++)
        {
                omegamats(i,0) = as<mat>(Omega[i]);
                omegapluses(i,0) = as<rowvec>(omegaplus[i]);
        }
        /* set sizes for sampling and return objects */
        alphacar = randu<rowvec>(numcar); /* sampling object */
        taucar = 0.1*ones<rowvec>(numcar); /* sampling object */
        Alphacar.set_size(nkeep,numcar); /*return object */
        Taucar.set_size(nkeep,numcar); /*return object */
	// Initialize log probability for Metropolis sampling alphacar
	logf.set_size(numcar);
        qmats.set_size(numcar,1);
    }
    
    if(numind > 0)
    {
        // fill in each field member with a vector of tau values for that type */  
        taus.set_size(numind,1);  /* sampling object */
        Tauind.set_size(numind,1);
        for(i = 0; i < numind; i++)
        {
            int dosei = numt(posind(i)); /* number of treatments for type i */
            rowvec tauind = 0.1*ones<rowvec>(dosei);
            taus(i,0) = tauind; 
            Tauind(i,0).set_size(nkeep,dosei); /* return object */
        }
    }
    // Put together P, treatment covariance matrix for delta
    mat Ptr(ntr,ntr); 
    buildP(pmatsmvn, taus, alphacar, taucar, omegamats, omegapluses,
           typet, Ptr); 
    
    // Build covariance matrix, Lambda, for the number of random effects
    mat lambda = 0.1*eye(nr,nr); /* sampling object */
    rowvec lambdavec(1,(nr*nr)); /* vectorize lambda to return in Lambda */
    mat Lambda(nkeep,(nr*nr)); /* return object */
    
    // Finally compose the covariance matrix for the nr*ntr x 1, delta(i)
    mat Pdelt = kron(lambda,Ptr); /* (nr*ntr) x (nr*ntr) */
    
    // Model Error
    colvec Taue(nkeep); /* model error */
    double taue = 1; /* sampling object */
    
    // DP return objects
    icolvec numM(nkeep);
    imat S(nkeep,np);
    field<icolvec> Num(nkeep,1);
    field<ucolvec> bigS;
    field<mat> optPartition(3,1); /* Hold empirical probability of co-clustering matrix, index of L-sq clustering scores objects, and cluster identifiers per iteration */

    // Initialize parameter values   
    /* cluster capture variables */
    int reps = np/M; /* truncato resulto */
    IntegerVector clust = seq_len(M); /* 1:M */
    IntegerVector s_ = rep_each(clust,reps);  /* 111222...MMM */
    icolvec s	= as<icolvec>(s_);
    icolvec srest((np-M*reps)); srest.ones(); srest *= M;
    s	= join_cols(s,srest); s -= 1; /* use s as position selector */
    s 	= shuffle(s,0);
    icolvec num(M);
    for(int i = 0; i < M; i++)
    {
	uvec pos	= find(s == i);
	num(i)          = pos.n_elem;
    }  
    // set locations for matrix of vectorized delta parameters, each nr*ntr
    mat    dstarmat(M,(nr*ntr)); rowvec mu(nr*ntr); mu.zeros(); rmvnrnd(Pdelt,dstarmat,mu);
    mat    dmat(np,(nr*ntr)); dmat.zeros();
    mat    thetamat(np,nr); thetamat.zeros();
    colvec zd(nc); zd.zeros();
    double conc = 1; /* DP concentration parameter */
    /* remaining parameters */
    colvec beta = randn<colvec>(nf)*sqrt(10);
    colvec resid(nc); resid.zeros(); /* regression residual */
    double alpha = rnorm( 1, 0, sqrt(10) )[0];
    /* fit assessment - related measures */
    double deviance = 0;  ucolvec ordscore; ordscore.zeros();
    mat phat(np,np); phat.zeros(); /* empirical co-clustering probability matrix */
    rowvec devmarg(nc); colvec devres(4); devres.zeros();
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    
    /* capture samples to return */
    int oo = 0, kk;

    // set hyperparameter values
    double a6, b6;
    a6 = ac; b6 = bc; /* conc */
    
    // Conduct posterior samples and store results
    int k,l, m;
    for(k = 0; k < niter; k++)
    {
	//if( (k % 1000) == 0 ) Rcout << "Interation: " << k << endl;
        clusterddpstep(zsplit, xsplit, ysplit, ytilsplit, doseperson, Pdelt, 
		       zbydose, dosequad, beta, alpha, taue, dstarmat, s, num, M, conc);
        concstep(conc, M, np, a6, b6);
        dstarstep(ytilsplit, zbydose, dosequad, taue, Pdelt, doseperson, s, num, dstarmat,
            zd, dmat, thetamat, nr);
        betaddpstep(xmat, xtx, y, beta, alpha, taue, zd, nf);
        alphaddpstep(xmat, beta, zd, y, resid, alpha, taue, nc);
        precisionddpstep(pmatsmvn, taus, lambda, alphacar, taucar, Ptr, Pdelt, taue,
                dstarmat, omegamats, omegapluses, resid, typet, numt, 
                logf, qmats, nc);
        if(k >= nburn)
        {
            deviance =  dev(resid, taue); /* scalar double deviance */
            dmarg(resid, taue, devmarg); /* 1 x nc vector of densities*/
            kk = k - nburn;
            if( kk == ((oo+1)*nthin - 1) )
            {
                Deviance(oo) = deviance;
                Devmarg.row(oo)  = devmarg;
                Beta.row(oo) = trans(beta);
                Alpha(oo) = alpha;
                Conc(oo) = conc;
                // return np*nr theta and b samples
                for( l = 0; l < nr; l++)
                {
                        thetarow.cols( np*l,(np*(l+1)-1) ) = trans( thetamat.col(l) );
                        /* rows of dmat fast index is treatment; slow index is order */
                        bint.cols( np*l,(np*(l+1)-1) ) = trans( dmat.col(ntr*l) ); 
                }
                Theta.row(oo) = thetarow;
		DoseEffects(nty,0).row(oo) = bint; /* random intercepts */
                // return np*(nr*dimm) dose effect samples for each treatment type, m
                // note: each "m" ties to treatment entry order in typet
                extractE(dmat, typet, numt, dmats, dosemats, dosevecs, nr);
                for(m = 0; m <nty; m++)
                {
                    dim = numt(m);
                    int dimp = nr*dim;
                    for(i = 0; i < np; i++)
                    {
                        doserow.set_size(np*dimp);
                        doserow.cols( dimp*i, (dimp*(i+1)-1) ) = trans( dosevecs(i,m) );
                    }
                    /* each row of DoseEffect(m,0) is an ntr_m treatment m dosage block for each of nr sets of effects */
                    /* for each of np subjects. So, for person i and effect order j, there are ntr_m dosages for treatment m*/
                    DoseEffects(m,0).row(oo) = doserow;
                }
                
                // return samples for precision parameters
                if(numcar > 0)
                {
                    Alphacar.row(oo) = alphacar;
                    Taucar.row(oo) = taucar;
                }
                if(numind > 0)
                {
                    for(i = 0; i < numind; i++)
                    {
                        Tauind(i,0).row(oo) = taus(i,0);
                    }
                }
                if(nummvn > 0)
                {
                    for(i = 0; i <nummvn; i++)
                    {
                        dim = pmatsmvn(i,0).n_rows;
			rowvec pmativec = reshape(pmatsmvn(i,0),1,(dim * dim)); /* by column */
                        Pmvn(i,0).row(oo) = pmativec; 
                    }
                }
                lambdavec = reshape(lambda,1,(nr*nr)); /* by column */
                Lambda.row(oo) = lambdavec;
                Taue(oo) = taue;
                Resid.row(oo) = trans(resid);
                numM(oo) = M;
                S.row(oo) = trans(s);
                Num(oo,0) = num;
                oo += 1; /* increment sample return counter */
            } /* end conditional statement on returning sample */
        } /* end if k > burnin, record results for return */
    } /* end MCMC loop over k */
    CARprecision(0,0) = Alphacar;
    CARprecision(1,0) = Taucar;
    
    // compute least squares cluster
    lsqcluster(S, Num, ordscore, phat, bigS);
    optPartition(0,0) = phat;
    optPartition(1,0) = conv_to<mat>::from(ordscore);
    optPartition(2,0) = conv_to<mat>::from(S);
    // DIC
    dic3comp(Deviance, Devmarg, devres); /* devres = c(dic,dbar,dhat,pd) */
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  //Rcpp::Named("Devmarg") = Devmarg,
                                  Rcpp::Named("devres") = devres,
                                  //Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("C") = Conc,
                                  Rcpp::Named("DoseEffects") = DoseEffects,
                                  Rcpp::Named("Theta") = Theta,
                                  Rcpp::Named("CAR_Q") = CARprecision,
                                  Rcpp::Named("Tauind") = Tauind,
                                  Rcpp::Named("Pmvn") = Pmvn,
                                  Rcpp::Named("Lambda") = Lambda,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Residuals") = Resid,
                                  Rcpp::Named("M") = numM,
                                  //Rcpp::Named("S") = S,
                                  Rcpp::Named("Num") = Num,
                                  Rcpp::Named("optPartition") = optPartition,
                                  Rcpp::Named("bigSmin") = bigS
				  );

END_RCPP
} /* end MCMC function returning SEXP */

   SEXP buildQ(const mat& omegamat, const rowvec& omegaplus, double alpha,
                double tau, mat& qmat, int& dim) 
    {
        BEGIN_RCPP
        // Function to build CAR precision matrix, qmat =  tau(D - alpha*Omega)
        dim = omegamat.n_rows;
        qmat.set_size(dim,dim); 
        qmat.zeros();
        qmat.diag() = omegaplus;
        qmat += -alpha*omegamat;
        qmat *= tau;
       END_RCPP
    } /* end function to build sparse nr x nr*ntr dosage matrix for each client*/
   
   SEXP extractE(const mat& dstarmat, const icolvec& typet, const icolvec& numt,
                field<mat>& dstarmats, field<mat>& cstarmats, 
                field<colvec>& cstarvecs, int nr) 
    {
        BEGIN_RCPP
        int M = dstarmat.n_rows; /* includes intercept parameters, mu (nr x 1) */
        int ntr = dstarmat.n_cols / nr;
        int nty = typet.n_elem;
        int m, k, dim;
        // create dstarmats
        rowvec dstarslice(nr*ntr); mat dslicemat(nr,ntr);
        // memo: dstarmats(M,0); cstarmats(M,nty); cstarvecs(M,nty);
        for(m = 0; m < M; m++)
        {
            dstarslice = dstarmat.row(m); /* 1 x ntr*nr */
            dslicemat = trans( reshape(dstarslice,ntr,nr) );
            dstarmats(m,0) = dslicemat;  
        }

        int startcol = 1, endcol; /*col 0 is the intercept */
        for(k = 0; k < nty; k++)
        {
            dim = numt(k);
            endcol = startcol + dim - 1;
            for(m = 0; m < M; m++)
            {
                cstarmats(m,k) = dstarmats(m,0).cols(startcol,endcol);
                //cstarvecs(m,k) = reshape(cstarmats(m,k),nr*dim,1,1); /* by row to get set of dosages */
                cstarvecs(m,k) = reshape(cstarmats(m,k).t(),nr*dim,1);
            }
            startcol += dim;   
        }
       END_RCPP
    } /* end function to build sparse nr x nr*ntr dosage matrix for each client*/

   SEXP buildP(const field<mat>& pmatsmvn, const field<rowvec>& taus, 
           const rowvec& alphacar, const rowvec& taucar, 
           const field<mat>& omegamats, const field<rowvec>& omegapluses,
           const icolvec& typet, mat& Ptr)
    {
        BEGIN_RCPP
        // Function to build covariance matrix for Treatment dosages
        // Composed of block diagonal structure with members a particular
        // covariance type determined by 'typet'.
        int m; /* loop index over treatment types */
        int sel = 0, dim = 0; /* dimension of local matrix for typet(m) */
        Ptr.zeros(); /* set elements to zero to build block diagonal */
        Ptr(0,0) = 1; /* intercept for Ptr, the treatments covariance matrix */
        int countmvn = 0, countind = 0, countcar = 0;
        int nty = typet.n_elem;
        mat qmat, pmat;
        int startP = 1, endP; /* first position is 1 for intercept */
        for(m = 0; m < nty; m++)
        {
            if(typet(m) == 1) /* CAR */
            {
               sel = countcar;
               buildQ(omegamats(sel,0),omegapluses(sel,0),alphacar(sel),
                       taucar(sel),qmat,dim);
               endP = startP + (dim - 1);
               Ptr.submat(span(startP,endP),span(startP,endP)) = qmat;
               countcar += 1;
            }else{
                if(typet(m) == 2) /* MVN */
                {
                    sel = countmvn;
                    dim = pmatsmvn(sel,0).n_rows;
                    endP = startP + (dim - 1);
                    pmat = pmatsmvn(sel,0);
                    Ptr.submat(span(startP,endP),span(startP,endP)) = pmat; 
                    countmvn += 1;
                }else{ /* typet(m) = 3 (IND) */
                    sel = countind;
                    dim = taus(sel,0).n_elem;
                    endP = startP + (dim - 1);
                    pmat.set_size(dim,dim); pmat.zeros();
                    pmat.diag() = taus(sel,0);
                    Ptr.submat(span(startP,endP),span(startP,endP)) = pmat; 
                    countind += 1;
                } /* end conditional loop on whether typet(m) is 2 or 3 */
            } /* end conditional statement on typet(m) */ 
            startP += dim; /* reset starting position for writing next entry in Ptr */
                    
        } /* end loop over m, the number the treatment types */
        END_RCPP
    } /* end function to compose treatments covariance matrix, Ptr */
   
    SEXP buildD(field<mat>& doseperson, const mat& dosemat, int nr) 
    {
        BEGIN_RCPP
        // Function to build the D_i, the nr x (nr*ntr), dosage design
        // matrix for each person which is back multiplied by the (nr*ntr) vectorized
        // dosage effect parameters, dstarmat(m) or dmat(i) and front multiplied
        // by random effects covariance matrix for subject i (or cluster m), zi. 
        int np = doseperson.n_rows;
        int ntr = dosemat.n_cols; /* dosemat is of dim = np x ntr */
        int i, k; mat dmat(nr,(nr*ntr)); 
        for(i = 0; i < np; i++)
        {
            dmat.zeros();
            int startcol = 0, endcol;
            for(k = 0; k < nr; k++)
            {
                endcol = startcol + (ntr - 1);
                dmat( k,span(startcol,endcol) ) = dosemat.row(i); /* person doses */
                startcol += ntr;
            }
            doseperson(i,0) = dmat;    
        }
       END_RCPP
    } /* end function to build sparse nr x nr*ntr dosage matrix for each client*/

    SEXP zdcomp(const field<mat>& doseperson, const mat& dmat, const mat& zmat,
            mat& thetamat, icolvec& persons, colvec& zd)
    {
        BEGIN_RCPP
        // Function to compute the product of z_i*D_i*delta_i for the client
        // random effects, producing an nc x 1 vector.
        int np = doseperson.n_rows;
        // compute number of repeated measures per person
        // positions is np x 1 with each cell the number of measures
        icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }
        // compute zd, nc x 1
        int startrow = 0;
        int j, nj, endrow; mat zj;
        for(j=0; j < np; j++)
        {
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            thetamat.row(j) = trans( doseperson(j,0)*trans( dmat.row(j) ) );
            zd.rows(startrow,endrow) = zj*trans( thetamat.row(j) );
            //re-set start positions
            startrow += nj;
        } /* end loop j for constructing thetamat and zd*/
        END_RCPP
    } /* end function to compute zd*/
    
    SEXP logp(const mat& omegamat, const rowvec& omegaplus,
            const field<colvec>& cstarvecs, const mat& lambda, double tau, 
            double alpha, int treat, double& logf, mat& qmat) 
    {
        BEGIN_RCPP
        // log probability calcuation for alpha parameter under proper CAR prior
        // 'treat' selects the focus treatment in typet
        int M = cstarvecs.n_rows;
        int nr = lambda.n_rows;
        int m, dim;
        buildQ(omegamat,omegaplus,alpha,tau,qmat,dim);
        qmat /= tau; /* CAR precision matrix */
        mat alphaprec = symmatl(kron(lambda,omegamat));
               
        double alphaterm = 0;
        for(m = 0; m < M; m++)
        {
               colvec bstarm = cstarvecs(m,treat);
               alphaterm += as_scalar( trans(bstarm)*alphaprec*bstarm );
                          
        }
        alphaterm *= 0.5*alpha*tau;
        // Probability of move for alpha
        double sign;
        double val;
        log_det(val, sign, qmat);
        logf = 0.5*nr*M*(val) + alphaterm;
        END_RCPP
    } /* end function to build sparse nr x nr*ntr dosage matrix for each client*/
    
    SEXP initlogp(const field<mat>& omegamats, const field<rowvec>& omegapluses,
            const field<colvec>& cstarvecs, const mat& lambda, const rowvec& taucar, 
            const rowvec& alphacar, const icolvec& typet, colvec& logfstart, 
            field<mat>& qmats) 
    {
        BEGIN_RCPP
        // Initial log prob calculation for vector of alpha parameters under proper CAR priors
        int nty = typet.n_elem;
        int k; /* loop over treatment type */
        int sel; /* focus treatment type using CAR covariance */
        int countcar = 0;
        double logf;
        mat omegamat, qmat;
        rowvec omegaplus;
        for(k = 0; k < nty; k++)
        {
            if(typet(k) == 1) /* CAR treatment */
            {
                sel = countcar;
                omegamat = omegamats(sel,0);
                omegaplus = omegapluses(sel,0);
                logp(omegamat, omegaplus, cstarvecs, lambda, taucar(sel), alphacar(sel), 
                        k, logf, qmat); 
                logfstart(sel) = logf;
                qmats(sel,0) = qmat;
                countcar += 1;
            }
                
        }       
        END_RCPP
    } /* end function to build sparse nr x nr*ntr dosage matrix for each client*/
    
    SEXP precisionddpstep(field<mat>& pmatsmvn, field<rowvec>& taus, mat& lambda,
           rowvec& alphacar, rowvec& taucar, mat& Ptr, mat& Pdelt, double& taue,
           const mat& dstarmat, const field<mat>& omegamats, 
           const field<rowvec>& omegapluses, const colvec& resid, 
           const icolvec& typet, const icolvec& numt, 
           colvec& logf, field<mat>& qmats, int nc)
    {
        BEGIN_RCPP
        // Function to build covariance matrix for Treatment dosages
        // Composed of block diagonal structure with members a particular
        // covariance type determined by 'typet'.
        int m; /* loop index over clusters */
        int k; /* loop index over treatment types */
        int d; /* loop over dosages within a treatment type */
        int sel = 0, dim = 0; /* dimension of local matrix for typet(k) */
        int countmvn = 0, countind = 0, countcar = 0;
        int nty = typet.n_elem;
        int M = dstarmat.n_rows;
        int nr = lambda.n_rows;
        field<mat> dstarmats(M,1); /* nr x ntr effects matrix for each of M clusters */
        field<mat> cstarmats(M,nty); /* nr x dim_k effects matrix for trt type k for each of M clusters */
        field<colvec> cstarvecs(M,nty); /* nr*dim_k x 1 effects vector for trt type k for each of M clusters */
        extractE(dstarmat, typet, numt, dstarmats, cstarmats, cstarvecs, nr);
        
        // sample treatment precision parameters
        for(k = 0; k < nty; k++)
        {
            dim = numt(k);
            if(typet(k) == 1) /* CAR */
            {
               sel = countcar;
               // Metropolis move for alpha
               /* draw candidate for alpha */
               double alpha = runif(1,0.6,1)[0];
               mat qmatnew, qmatold; /* CAR precision matrix */
               mat omegamat = omegamats(sel,0);
               rowvec omegaplus = omegapluses(sel,0);
               double tau = taucar(sel);
               double logfnew, logfold;
               // compute log probabilities for current and proposed values for alpha
	       logp(omegamat, omegaplus, cstarvecs, lambda, tau, alphacar(sel), k, logfold, qmatold);
               logp(omegamat, omegaplus, cstarvecs, lambda, tau, alpha, k, logfnew, qmatnew); 
               // compute probability of move
               double logp = logfnew - logfold;
               double thresh = runif(1,0,1)[0];
               if(log(thresh) < logp)
               {
                   alphacar(sel) = alpha;
                   logf(sel) = logfnew;
                   qmats(sel) = qmatnew;
               }else{
		   logf(sel) = logfold;
		   qmats(sel) = qmatold;
               }        
               // Gibbs sample for tau
               double tauterm = 0;
               mat tauprec   = symmatl(kron(lambda,qmats(sel)));
               for(m = 0; m < M; m++)
               {
                   colvec bstarm = cstarvecs(m,k);
                   tauterm  += as_scalar( trans(bstarm)*tauprec*bstarm );   
               }
               tauterm *= 0.5;
               tauterm += 1;
               double astar = 0.5*nr*dim*M + 1;
               taucar(sel) = rgamma(1, astar, (1/tauterm))[0];
               // increment counter to next CAR treatment (if one exists)
               countcar += 1;
            }else{
                if(typet(k) == 2) /* MVN */
                {
                    sel = countmvn;
                    double nu = dim + 1;
                    double nupost = nu + nr*M;
                    mat pmatpost = 0.1*eye(dim,dim);
                    for(m = 0; m < M; m++)
                    {
                        pmatpost += trans(cstarmats(m,k))*lambda*cstarmats(m,k);
                    }
                    mat vmatpost = inv(pmatpost);
                    mat pnew(dim,dim);
                    wishrnd(pnew, vmatpost, nupost);
                    pmatsmvn(sel,0) = pnew;
                    countmvn += 1;
                }else{ /* typet(k) = 3 (IND) */
                    sel = countind;
                    rowvec taunew(dim);             
                    double bstar, astar;
                    astar = 0.5*nr*M + 1;
                    for(d = 0; d < dim; d++)
                    {
                        bstar = 0;
                        for(m = 0; m < M; m++)
                        {
                           colvec cstarmd = cstarmats(m,k).col(d);
                           bstar += as_scalar( trans(cstarmd)*lambda*cstarmd );
                        }
                        bstar *= 0.5;
                        bstar += 1;
                        taunew(d) = rgamma(1, astar, (1/bstar))[0];
                    }
                    taus(sel,0) = taunew;
                    countind += 1;
                } /* end conditional loop on whether typet(k) is 2 or 3 */
            } /* end conditional statement on typet(k) */                    
        } /* end loop over k, the number the treatment types */
        
        // sample random effect precision matrix, lambda
        buildP(pmatsmvn, taus, alphacar, taucar, omegamats, omegapluses,
           typet, Ptr); /* compose new Ptr from sampled treatment precisions */
        double nu = nr + 1;
        double nupost = nu + M*sum(numt);
        mat lambdapost = 0.1*eye(nr,nr);
        for(m = 0; m < M; m++)
        {
              lambdapost += dstarmats(m)*Ptr*trans(dstarmats(m));
        }
        mat vlampost = inv(lambdapost);
        wishrnd(lambda, vlampost, nupost);
        // compose new Pdelt for sampling dstarmat
        Pdelt = kron(lambda,Ptr);
        
        // sample error precison, taue
        double a = 0.5*double(nc) + 1;
        double b = 0.5*dot( resid, resid ) + 1;
        taue = rgamma(1, a, (1/b))[0];
        
        END_RCPP
    } /* end function to sample all model precision parameters */
    
    SEXP clusterddpstep(const field<mat>& zsplit, const field<mat>& xsplit, const field<colvec>& ysplit,
            field<colvec>& ytilsplit, const field<mat>& doseperson, const mat& Pdelt, 
	    const field<mat>& zbydose, const field<mat>& dosequad,
            const colvec& beta, double alpha, double taue, 
            mat& dstarmat, icolvec& s, icolvec& num, int& M, 
            double& conc)
    {
        BEGIN_RCPP
        int np = doseperson.n_rows;
        int nr = zsplit(0,0).n_cols;
        int ntr = dstarmat.n_cols / nr;
        // compute number of cases for each person and store results in
        
        // sample cluster assignments, s(1), ..., s(np)
        // fixing all other parameters, including dstarmat
        mat phid(nr*ntr,nr*ntr);
        colvec cj, cstarj, ytildej, dstarj(nr*ntr), ed(nr*ntr), hd(nr*ntr);
        double logq0, q0, sweights;
        int j, l, m;
        for(j = 0; j < np; j++)
        {
            /* subtracting cj eliminates need to t-fer xj to sample bstar */
	    cj			= alpha + xsplit(j,0)*beta;
	    ytildej 		= ysplit(j,0) - cj;
            ytilsplit(j,0) 	= ytildej;

            // sample cluster assignment indicators, s(1),...,s(np)
            if(num(s(j)) == 1) /* remove singleton cluster */
            {
                dstarmat.shed_row(s(j));
                num.shed_row(s(j));
                M -= 1; /* decrement cluster count */

                //decrement cluster tracking values by 1 for tossed cluster
                uvec ord = sort_index(s); /* of size np, same as s */
                int k; int thresh = np;
                for (k = 0; k < np; k++)
                {
                    if( s(ord(k)) > s(j) )
                    {
                        thresh = k;
                        break;
                    }
                } /* end loop k to find thresh*/
                // decrement s > s(j)
                for (m = thresh; m < np; m++)
                {
                    s(ord(m)) -= 1;
                }
            } /* end cluster accounting adjustment for singleton cluster */
            else
            {
                num(s(j)) -= 1;
            } /* decrement non-singleton cluster count by one */

            // construct normalization constant, q0, to sample s(j)
            // compute posterior mean and variance to compute density for bj
            // build loqq0 and exponentiate
            logq0 = loglike(ytildej,taue); /* likelihood */
            colvec zro(nr*ntr); zro.zeros(); /* add prior with zero mean */
            logq0 += logdens(zro,Pdelt);
            // build posterior density
            ed = taue*zbydose(j,0).t()*ytildej;
            phid = taue*dosequad(j,0) + Pdelt; /*nr*ntr x nr*ntr */
	    hd = solve(phid,ed);
            logq0 -= logdens(hd,phid); /* dmvn(hd,phibd^-1) */
            q0 = exp(logq0);

            // construct posterior sampling weights to sample s(j)
            colvec weights(M+1); weights.zeros();
            for (l = 0; l < M; l++) /* cycle through all clusters for s(j) */
            {
                s(j) 	= l; /* will compute likelihoods for every cluster */
                dstarj 	= trans( dstarmat.row(s(j)) ); /* nr*ntr */
                cstarj 	= cj + zbydose(j,0)*dstarj;
                ytildej = ysplit(j,0) - cstarj;

                weights(l) = exp(loglike(ytildej,taue));
                weights(l) *= num(s(j))/(double(np) - 1 + conc);
            }
            weights(M) = conc/(double(np)-1+conc)*q0;

            // normalize weights
            sweights = sum(weights);
            if(sweights == 0)
            {
                weights.ones(); weights *= 1/(double(M)+1);
            }
            else
            {
                weights = weights / sweights;
            }

            // conduct discrete posterior draw for s(j)
            unsigned long MplusOne = M + 1;
            s(j) = rdrawone(weights, MplusOne);

            // if new cluster chosen, generate new location
            if(s(j) == M)
            {
                dstarj.zeros();
                rmvnchol(zbydose(j,0), dosequad(j,0), Pdelt, ysplit(j,0), cj, dstarj, nr*ntr, taue);
                dstarmat.insert_rows(M,1);
                num.insert_rows(M,1);
                dstarmat.row(M) = trans(dstarj);
                num(M) = 1;
                M = MplusOne;
            }
            else
            {
                num(s(j)) += 1;
            }
            
        } /* end loop j for cluster assignment */
        END_RCPP
    } /* end function dstep for cluster assignments, s. */

    
    SEXP dstarstep(const field<colvec>& ytilsplit,
	    const field<mat>& zbydose, const field<mat>& dosequad,
            double taue, const mat& Pdelt, const field<mat>& doseperson, 
            const icolvec& s, const icolvec& num, mat& dstarmat,
            mat& zd, mat& dmat, mat& thetamat, int nr)
    {
        BEGIN_RCPP
        // sample dmat, independently, by person
        int M = dstarmat.n_rows;
        int np = doseperson.n_rows;
        int ntr = dstarmat.n_cols / nr;
        int j, k, l, m, numj;
        colvec dstarj(nr*ntr);
        
        // collect subsets of cluster indices, i: s(i) = j, j = 1, ..., M

        uvec ord = sort_index(s); /* of length np > = M*/
        uvec thresh(M); thresh(0) = 0;
        int count = 1;
        //find cluster index start (change) points
        for (m = 0; m < (np-1); m++)
        {
            if( s(ord(m+1)) > s(ord(m)) )
            {
                thresh(count) = m+1;
                count += 1;
            }
        }
        // sample dstarmat for each cluster, j
        int rowchoose; 
        for(j=0; j < M; j++) /* j is now the cluster label */
        {
            numj = num(j);
            field<mat> zjwedge(numj,1);
	    field<mat> zquadjwedge(numj,1);
            field<colvec> yjtilwedge(numj,1);
            // collect data observations - z,ytilde - for each cluster
            for(k = 0; k < numj; k++)
            {
                /* extract by-person matrix views from data*/
                rowchoose 		= ord(thresh(j) + k);
		zjwedge(k,0) 		= zbydose(rowchoose,0);
		zquadjwedge(k,0)	= dosequad(rowchoose,0);
                yjtilwedge.row(k) 	= ytilsplit.row(rowchoose);
            }
            /* sample posterior of nj block of dstarmat*/
            dstarj.zeros();
            rmvncholclust(zjwedge, zquadjwedge, Pdelt, yjtilwedge, dstarj, taue);
            dstarmat.row(j) = trans(dstarj);

        } /* end loop j for sampling dstarj */

        // construct dmat for function return
         int startrow = 0; int endrow, nl; rowvec dlr(nr*ntr);
         for(l = 0; l < np; l++) /* now we're iterating by person */
            {
                /* dmat */
                dlr = dstarmat.row(s(l));
                dmat.row(l) = dlr;
                /* zd */
                nl = zbydose(l,0).n_rows;
                endrow = startrow + nl - 1;
                zd.rows(startrow,endrow) = zbydose(l,0)*trans(dlr);
                startrow += nl;
                /* thetamat (np x nr) */
                thetamat.row(l) = trans( doseperson(l,0)*trans(dlr) );
            }

        END_RCPP
    } /* end function dstarstep for sampling dmat and zd */

    SEXP betaddpstep(const mat& xmat, const mat& xtx, const colvec& y,
                colvec& beta, double alpha, double taue,
                const colvec& zd, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zd; /* nc x 1 */
        rmvnlschol(xmat, xtx, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */


    SEXP alphaddpstep(const mat& xmat, const colvec& beta,
                const colvec& zd, const colvec& y,
                colvec& resid, double& alpha, double taue,
                long nc)
    {
        BEGIN_RCPP
        colvec c = xmat*beta + zd;
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc);
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */

    SEXP CovUnitSplitDDP(const mat& zmat, const mat& xmat, const field<mat>& doseperson, const colvec& y, field<mat>& zsplit, 
			field<mat>& xsplit, field<colvec>& ysplit, field<mat>& dosequad, field<mat>& zbydose, icolvec& persons)
   {
	BEGIN_RCPP
	// initialize objects
	int np = zsplit.n_rows; 
	int startrow, endrow, nj, j;
	mat zj, xj; colvec yj;

	// compute number of repeated measures for each unit.
 	// assumes repeated units listed in nested / block fashion under unit
	icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }

	startrow	= 0;
    	for(j = 0; j < np; j++)
        {
            // store by-person data cuts for later sampling of clust locs
            /* extract by-person matrix views from data*/
            nj 			= positions(j);
            endrow 		= startrow + nj - 1;
            zj 			= zmat.rows(startrow,endrow);
            xj 			= xmat.rows(startrow,endrow);
	    yj 			= y.rows(startrow,endrow);
	    zbydose(j,0)	= zj * doseperson(j,0);
	    dosequad(j,0)	= zbydose(j,0).t() * zbydose(j,0);
            zsplit(j,0) 	= zj; /* np rows */
	    xsplit(j,0) 	= xj;	    
	    ysplit(j,0) 	= yj;
	    startrow		+= nj;
	}
	END_RCPP
    } /* end function splitting design objects by clustering unit */
    
