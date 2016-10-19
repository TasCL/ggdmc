//   Copyright (C) <2016>  <Yi-Shin Lin>
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; version 2
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include <regex>
#include <ggdmc.hpp>
#include <RcppArmadillo.h>

// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP; eg OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(cpp11)]]
bool checkDDM(std::vector<double> pVec) {
  /* pVec entering checkDDM is in Voss's form before appending RT. pVec does
  not distinguish different factor levels eg v.f1, v.f2. It is ordered as
  ----------Name sequence correction----------
  * 0  1  2  3  4   5   6  7    8   9
  * a  v  z  d  sz  sv  t0 st0  RT  precision (dmc uses names)
  * Combine Voss, Andrew, Matthrew checks all in one */
  //double aLimit = 2.0 ;  // Use Andrew's limits; Voss=10.0; Verdonck=1.0
  //double vLimit = 5.0 ;  // Voss=50.0
  //double zr = pVec[2] / pVec[0] ;
  double szLower  = pVec[2] - 0.5*pVec[4] ;
  double szUpper  = pVec[2] + 0.5*pVec[4] ;
  double szLimit = 1.999999*std::min(pVec[2], 1-pVec[2]) ;

  // In order to produce exact same profile, do less checks
  // This will make sampler less efficient, because making more failure guesses
  // 1. a cannot be less than z, or less than (or equal to) 0 or greater than
  // its limit
  bool c1 = pVec[0] < 0 ;
  // 2. If v's absolute value greater than vLimit (Voss's 50; Andrew's 5)
  // bool c2 = std::abs(pVec[1]) > vLimit ;
  // 3. Check zr. If zLower less than 1e-6 or zUpper greater than 1
  bool c3 = (szLower < 1e-6 || szUpper > 0.999999) ;  // 0 ~ 1
  // 4. If t0 - abs(d)/2 - st0/2 is less than 0
  bool c4 = (( pVec[6] - 0.5*std::abs(pVec[3]) - 0.5*pVec[7] ) < 0) ;
  // 5. If t0 is almost 0
  bool c5 = pVec[6] < 1e-6 ;
  // 6. If szr greater than 1 or less than 0 or greater than 2 x min(z)
  bool c6 = (pVec[4]  < 0) || pVec[4] > szLimit ;
  // 7. If sv greater than 2 or less than 0
  bool c7 = pVec[5] > 2 || pVec[5] < 0;
  // 8. If st0 less than 0
  bool c8 = pVec[7] < 0 ;
  // 9. If t0 - st0/2 less than 0
  bool c9 = (pVec[6] - pVec[7]/2) < 0 ;
  bool out = c1 || c3 || c4 || c5 || c6 || c7 || c8 || c9 ;

  return out;
}

//' Map a parameter vector to an accumulator matrix
//'
//' An Rcpp function to get accumulator x parameter matrix for a design cell. For
//' two accumulator model, like drift-diffusion model, the accumulators are r1
//' and r2. For mulitple-accumulator model, like LBA model, the accumulators
//' are r1, r2, r3, etc. But because LBA probability density has yet
//' implemented, the function now only works for DDM. The part calling
//' \code{transform.dmc} is replaced by ordering the parameter sequence inside
//' \code{density.cpp}. The \code{n1order} switch is still kept, until the R's
//' \code{simulate.dmc} is upgraded to C++ version. This function is akin to
//' DMC's \code{p.df.dmc}. This is a Rcpp function.
//'
//' @param pVec a user supplied parameter vector or a sampler supplied theta/phi
//' vector.
//' @param cell a string indicating the name of a experimental condition
//' @param model an EAM model specification, created usually by model.dmc. Read
//' as a NumericVector, so attributes are kept.
//' @param n1order a boolean value indicating if swap accumulator row to use
//' "n1" order. This is for LBA model
//' @return an accumulator x EAM parameter matrix (without factor level).
//' @export
//' @examples
//' m1 <- model.dmc(
//'     p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="F",t0="1",st0="1"),
//'     match.map = list(M=list(s1="r1",s2="r2")),
//'     factors=list(S=c("s1","s2"), F=c("f1","f2")),
//'       constants = c(st0=0,d=0),
//'       responses = c("r1","r2"),
//'       type = "rd")
//' # Parameter vector names are: ( see attr(,"p.vector") )
//' # [1] "a"     "v.f1"  "v.f2"  "z"     "sz"    "sv.f1" "sv.f2" "t0"
//' #
//' # Constants are (see attr(,"constants") ):
//' #   st0   d
//' #     0   0
//' #
//' # Model type = rd
//'
//' attr(m1,"p.vector")
//' # a  v.f1  v.f2     z    sz sv.f1 sv.f2    t0
//' # NA    NA    NA    NA    NA    NA    NA    NA
//' pVec <- c(1,1,1,1, 1,1,1,1)
//' names(pVec) <- names(attr(m1, "p.vector"))
//'
//' accMat1 <- getAccumulatorMatrix(pVec, "s1.f1.r2", m1, FALSE)
//' ##         col_names
//' ##row_names a v z d sz sv t0 st0
//' ##       r1 1 1 0 0  1  1  1   0
//' ##       r2 1 1 0 0  1  1  1   0
//' accMat2 <- getAccumulatorMatrix(pVec, "s1.f1.r2", m1, TRUE)
// [[Rcpp::export]]
Rcpp::NumericMatrix getAccumulatorMatrix(Rcpp::NumericVector pVec,
  std::string cell, Rcpp::NumericVector model, bool n1order=true) {
  // Get all parameters NA/0 vector. These include eg v.f1, v.f2 as well as
  // d, st0 etc, the parameters set by default.
  Rcpp::NumericVector allParameters   = model.attr("all.par") ; // NA, 0
  std::vector<std::string> allParName = allParameters.names() ;
  // Get the parameters set directly by the user, ie p.vector. They are ega,
  // v.f1, & v.f2, but not d & st0. The names in the pVecNA should be identical
  // as the p.vector set by prior.p.dmc
  Rcpp::NumericVector pVecNA      =  model.attr("p.vector") ;
  std::vector<std::string> pNames = pVecNA.names() ;

  // 1. Extract the 8 fastdm parameter names. This does not consider factor level
  // type.par.names gives identical names as par.names, but with a different
  // order
  // 2. type.par.names is used only in transform.dmc, but it is removed.
  // std::vector<std::string> fastdmNames2 = model.attr("type.par.names") ;
  std::vector<std::string> fastdmNames = model.attr("par.names") ; // 8


  /* Nidan: Build parameter matrix (in line with Andrew's par.mat) ---------
  * dim1 = (eg)
  [1] "s1.f1.r1" "s2.f1.r1" "s1.f2.r1" "s2.f2.r1" "s1.f1.r2"
  [6] "s2.f1.r2" "s1.f2.r2" "s2.f2.r2"
  * dim2 = "a"    "v.f1" "v.f2" "z"    "d"    "sz"   "sv"   "t0"   "st0"
  * dim3 = "r1" "r2"
  ------------------------------------------------------------------------- */
  Rcpp::List modelDimension     = model.attr("dimnames") ;
  std::vector<std::string> dim1 = modelDimension[0] ;  // row=cells
  std::vector<std::string> dim2 = modelDimension[1] ;  // col=allPars
  std::vector<std::string> dim3 = modelDimension[2] ;  // slice=accumulators
  // In case Andrew's n1.order uses a different ordering for the condition,
  // I will not reuse dim1.
  int rowLen = dim1.size() ;  // condition number, s1.r1 etc.
  int colLen = dim2.size() ;  // fast-dm parameter number
  int sliLen = dim3.size() ;  // accumulator number
  // Cast model to arma cube and name it as modelCube
  arma::cube modelCube(model.begin(), rowLen, colLen, sliLen, false) ;
  arma::mat modelMat ; // T/F conditions x all.pars
  Rcpp::NumericMatrix out(sliLen, fastdmNames.size()) ; // 2 x 8

  // Iterate through accumulators, eg r1, r2, (if have any), r3, etc.
  for(std::vector<std::string>::iterator itSli = dim3.begin();
    itSli != dim3.end(); ++itSli)
  {
    Rcpp::NumericVector accumulatorContainer(fastdmNames.size()) ; // 8
    int sliIdx = distance(dim3.begin(), itSli) ;

    // Iterate through conditions, eg s1.r1, s2.r1, s1.r2, s2.r2 etc.
    for (std::vector<std::string>::iterator itRow = dim1.begin();
      itRow != dim1.end(); ++itRow)
    {
      int rowIdx = distance(dim1.begin(), itRow) ;
      modelMat   = modelCube.slice(sliIdx) ; // ncells x nallpars
      // parameter boolean index, 'cause model is a boolean array
      // if this row is the indicated cell (condition;string)
      if (*itRow == cell)
      {
        int ii = 0 ;
        arma::vec cellBool = vectorise(modelMat.row(rowIdx)) ;
        // para1, para2, para3, (eg a, v.f1, v.f2, z, d, sz) ...
        for (std::vector<std::string>::iterator itCol = dim2.begin();
          itCol != dim2.end(); ++itCol) // 0-13
        {
          int colIdx = distance(dim2.begin(), itCol) ;
          // int parameterIndicator = cellBoolRow(colIdx) ;
          if (cellBool(colIdx)==1) // this presume rowSums = nAllPars
          {
            //--- Column name order uses "all.par" ----
            // a   v   z   d  sz  sv  t0 st0;
            accumulatorContainer[ii] = allParameters[colIdx] ;
            std::vector<std::string> pVecNames = pVec.names() ;
            for(int i=0; i<pVec.size(); i++)
            {
              bool isUsed = (pVecNames[i] == *itCol) ;
              if(isUsed == 1) { accumulatorContainer[ii] = pVec[i] ; }
            }
            ii++ ;
          }
        }
      }
    }
    out(sliIdx, Rcpp::_) = accumulatorContainer ;
  }

  /* -------------------------------------------------------------------------
  Third part - If model type equals to Ratcliff's DDM, flip parameter z
  ------------------------------------------------------------------------- */
  std::string modelType = model.attr("type") ; // Is it DDM, or other EAMs
  bool bool1 = (modelType == "rd") ;      // rd equals to DDM
  bool bool2 = false ;

  // is.r1 attribute is found in DDM model only
  if (bool1==true)
  {
    Rcpp::LogicalVector Response1Indicator = model.attr("is.r1") ;
    std::vector<std::string> ConditionNames = Response1Indicator.names() ;

    for (std::vector<std::string>::iterator it = ConditionNames.begin();
      it != ConditionNames.end(); ++it)
    {
      int idx = distance(ConditionNames.begin(), it) ;
      if (*it == cell ) { bool2 = Response1Indicator[idx]; } ;
    }
  }

  // 'cause the name of factor level have been removed, now I use fastdmNames,
  // instead of allParNames std::vector<std::string> fastdmNames
  if (bool1 && bool2)
  {
    for(std::vector<std::string>::iterator itSli = dim3.begin();
      itSli != dim3.end(); ++itSli)
    {
      int sliIdx = distance(dim3.begin(), itSli) ;
      for(std::vector<std::string>::iterator itfastdmNames=fastdmNames.begin();
        itfastdmNames != fastdmNames.end(); ++itfastdmNames)
      {
        int idx = distance(fastdmNames.begin(), itfastdmNames) ;
        // No need to consider "z.", 'cause now all in fastdm original names
        // substr still all right, because it takes "z" and "z." into account.
        std::string zDot = fastdmNames[idx].substr(0, 2) ;
        if (zDot == "z")
        {
          out(sliIdx, idx) = 1 - out(sliIdx, idx) ;
        }
      }
    }
  }

  /* -------------------------------------------------------------------------
  Fourth part - if n1order is true, swap accumulators
  ------------------------------------------------------------------------- */
  Rcpp::NumericMatrix res(sliLen, fastdmNames.size()) ;
  if (n1order) {
    // (row=condition) x (col=accumulator)
    Rcpp::NumericMatrix n1OrderMatrix = model.attr("n1.order") ;
    Rcpp::List n1OrderMatrixDimension = n1OrderMatrix.attr("dimnames") ;
    // n1Order's row names are a list of all condition names, eg s1.r1 ...
    std::vector<std::string> rowNames = n1OrderMatrixDimension[0] ;
    Rcpp::NumericVector targetCondContainer ;

    for (std::vector<std::string>::iterator itAllConditions = rowNames.begin();
      itAllConditions != rowNames.end(); ++itAllConditions)
    {
      if (*itAllConditions == cell)
      {
        int conditionIdx = distance(rowNames.begin(), itAllConditions) ;
        // Select the row matches to the cell name
        targetCondContainer = n1OrderMatrix(conditionIdx, Rcpp::_) ;
      }
    }

    for(std::vector<std::string>::iterator itSli = dim3.begin();
      itSli != dim3.end(); ++itSli)
    {
      int sliIdx = distance(dim3.begin(), itSli) ;
      res(sliIdx, Rcpp::_) = out(targetCondContainer[sliIdx]-1, Rcpp::_);
    }
  } else {
    res = out ;
  }

  //--- Column name order uses "all.par" ----
  // a   v   z   d  sz  sv  t0 st0;
  // This will be sorted later to fit my g_minus/g_plus function
  // Row name, by default, uses n1 order (ie set TRUE/true by default)
  // I put n1Order process firstly, so no need to negate bool n1order.
  Rcpp::List resNames = Rcpp::List::create(
    Rcpp::Named("row_names") = dim3,
    Rcpp::Named("col_names") = fastdmNames) ;

  res.attr("dimnames") = resNames ;
  return res ;
}


//' Compute Probability Density of Drift-Diffusion Model
//'
//' This function implements the equations in Voss, Rothermund, and Voss (2004).
//' These equations calculate Ratcliff's drift-diffusion model (RDM, 1978).
//' \code{ddmc} re-implements Voss, Rothermund, and Voss's (2004) equations A1
//' to A4 (page 1217) via Rcpp. This Rcpp function is akin to DMC's
//' \code{likelihood.dmc}. There is a \code{ddmc_parallel} with same arguments
//' using Open MPI to calculate RDM probability densities when data sets are
//' large (e.g., more than 5000 trials per condition) or when precision is set
//' high.
//'
//' @param x a data frame containing choices and RTs. This is akin to dnorm's x
//' argument.
//' @param pVec a parameter vector. For example,
//' p.vector <- c(a=1.25, v.f1=.20, v.f2=.33, z=.67, sz=.26, sv=.67, t0=.26)
//' @param precision a precision parameter. The larger, the higher precision to
//' estimate diffusion probability density. A general recommendation for the low
//' bound of precision is 2.5. If you need to use a precision higher than that,
//' you may want to consider using \code{ddmc_parallel}. A default value is set
//' at 2.5.
//' @param minLike a minimal log likelihood. If a estimated density is
//' less than minLike, the function returns minLike. A default value is set at
//' 1e-10.
//' @return a double vector with drift-diffusion probability density
//' @references Voss, A., Rothermund, K., & Voss, J. (2004).  Interpreting the
//' parameters of the diffusion model: An empirical validation.
//' \emph{Memory & Cognition}, \bold{32(7)}, 1206-1220. \cr\cr
//' Ratcliff, R. (1978). A theory of memory retrival. \emph{Psychological
//' Review}, \bold{85}, 238-255.
//' @export
//' @examples
//' ## Set up a basic RDM design
//' m1 <- model.dmc(
//'     p.map=list(a="1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
//'     constants=c(st0=0, d=0),
//'     match.map=list(M=list(s1="r1", s2="r2")),
//'     factors=list(S=c("s1", "s2")),
//'     responses=c("r1", "r2"),
//'     type="rd")
//'
//' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
//'
//' ## Set up a model-data instance
//' raw.data <- simulate(m1, nsim=1e2, p.vector=pVec)
//' mdi <- data.model.dmc(raw.data, m1)
//'
//' ## Print probability densities on the screen
//' ddmc(mdi, pVec)
//' ## [1]  0.674039828 0.194299671 0.377715949 1.056043402 0.427508806
//' ## [6]  2.876521367 0.289301835 2.701315781 0.419909813 1.890380685
//' ## [11] 2.736803321 1.278188278 0.607926999 0.024858130 1.878312335
//' ## [16] 0.826923114 0.189952016 1.758620294 1.696875882 2.191290074
//' ## [21] 2.666750382 ...
//'
//' mdi.large <- data.model.dmc(simulate(m1, nsim=5e3, p.vector=pVec), m1)
//'
//' ## d.large <- dplyr::tbl_df(mdi.large)
//' ## d.large
//' ## Source: local data frame [10,000 x 3]
//' ##
//' ##         S      R        RT
//' ##    (fctr) (fctr)     (dbl)
//' ## 1      s1     r2 0.2961630
//' ## 2      s1     r1 0.9745674
//' ## 3      s1     r2 0.6004945
//' ## 4      s1     r1 0.1937356
//' ## 5      s1     r1 0.3608844
//' ## 6      s1     r1 0.5343379
//' ## 7      s1     r1 0.2595707
//' ## 8      s1     r1 0.3247941
//' ## 9      s1     r1 0.4901325
//' ## 10     s1     r2 0.5011935
//' ## ..    ...    ...       ...
//'
//' system.time(den1 <- ddmc(mdi.large, pVec))
//' ##  user  system elapsed
//' ## 0.028   0.004   0.032
//'
//' system.time(den2 <- ddmc_parallel(mdi.large, pVec))
//' ##  user  system elapsed
//' ## 0.200   0.000   0.019
//' all(den1==den2)
//' ## [1] TRUE
// [[Rcpp::export]]
std::vector<double> ddmc(Rcpp::DataFrame x, Rcpp::NumericVector pVec,
  double precision=2.5, double minLike=1e-10)
{
  /* -------------------------------------------------------------------------
  First part - Create a boolean vector for choosing g_plus(1) or g_minus(0)
  ------------------------------------------------------------------------- */
  // Follow dnorm's sequence, I put quantile (named x) first (choice RTs) and
  // then the parameters (ie pVec); choice is achieved by ifelse "bound"
  // (see below). Others arguments eg fastdm's precision and GPU settings
  // come after (GPU setting has not added yet)
  Rcpp::NumericVector model = x.attr("model") ;         // For ddmc
  Rcpp::List matchMapOuter  = model.attr("match.map") ; // Handle two-layer list
  Rcpp::List matchMapInner  = matchMapOuter["M"] ;      // "M" as match.map's name
  Rcpp::List modelDimension = model.attr("dimnames") ;

  std::vector<std::string> stimulusName   = matchMapInner.names() ;
  std::vector<std::string> conditionNames = modelDimension[0] ; // row
  std::vector<std::string> parameterNames = modelDimension[1] ; // col; para
  std::vector<std::string> accumulatNames = modelDimension[2] ; // slice; r1, r2
  std::vector<int> bound (conditionNames.size());

  for (Rcpp::List::iterator it1 = matchMapInner.begin();
    it1 != matchMapInner.end(); ++it1)
  {
    int i = std::distance(matchMapInner.begin(), it1) ;
    std::string a = stimulusName[i] ;  // s1, s2
    std::string b = *it1 ;  // *it == matchMapInner[i]; r1, r2
    std::regex matchPattern (a + "(.*)" + b) ;
    for (std::vector<std::string>::iterator it2 = conditionNames.begin();
      it2 != conditionNames.end(); ++it2 )
    { // bound needs j
      int j = std::distance(conditionNames.begin(), it2) ;
      // Later accumulator test will overwrite previous, so test if only when
      // it is still not matched.
      if(bound[j] == 0)
      {
        bound[j] = (std::regex_match (*it2, matchPattern)) ? 1 : 0 ;
      }
    }
  }

  /* -------------------------------------------------------------------------
  Mawashigeri - min.like is set default in dmc.hpp
  ------------------------------------------------------------------------- */
  Rcpp::LogicalVector isempty = x.attr("cell.empty") ;
  Rcpp::List cellIdx = x.attr("cell.index") ; // s1.r1, s1.r2, s2.r1, s2.r2 etc.
  Rcpp::CharacterVector isemptyNames = isempty.names() ;
  arma::vec RT = x["RT"] ; // extract RTs from the input DataFrame (quantile)
  Rcpp::NumericMatrix accMat ; // forward declaration, avoid inside
  Rcpp::NumericVector accumulator ;  // loop.
  std::vector<double> out(RT.size()) ;  // ditto
  // prepare a container for pVec sending into g_minus/plus function, 'cause
  // my version of g_minus/plus require a redunctant element as RT
  std::vector<double> gVec(10) ;
  gVec.at(9) = precision ;

  arma::uvec selectedRTIdx ;  // 'cause arma's find function wants uvec type

  /* Let just use one vector, passing also RT and precision;
  * below is seq from dimnames(model)[[2]]
  * 0  1  2   3  4   5   6  7    8   9
  * a  v  zr  d  sz  sv  t0 st0  RT  precision (dmc uses names) */
  for (Rcpp::LogicalVector::iterator it3 = isempty.begin();
    it3 != isempty.end(); ++it3 )
  {
    // Why spend the effort to arrange the pVec for each accumulator and then
    // pick only the first accumulator?  Is it a bug? or just to accommodate rtdists?
    int k = std::distance(isempty.begin(), it3) ;
    accMat = getAccumulatorMatrix(pVec, conditionNames[k], model, true) ; // n1order==true-> LBA
    accumulator = accMat(0, Rcpp::_) ;

    // fill in the gVec with the 8 parameter values firstly
    int gVecLen = gVec.size() - 2 ;
    for(int ii = 0; ii < gVecLen; ii++) { gVec[ii] = accumulator[ii] ; }
    // Because cellIdx is a List, I can use name indexing
    // Select RTs in conditionNames[k], eg s1.r1, condition. In order to use
    // indexing, cast Rcpp's NumericVector to Arma's vec
    // > attr(data.model,"cell.empty")
    //   s1.r1 s2.r1 s1.r2 s2.r2
    //   FALSE FALSE FALSE FALSE
    Rcpp::NumericVector cppIdxRT = cellIdx[conditionNames[k]] ; // T, T, F, etc.
    arma::vec armaIdxRT(cppIdxRT.begin(), cppIdxRT.size(), false) ;
    selectedRTIdx = find(armaIdxRT == 1) ;  // the location index equal TRUE
    arma::vec selectedRT = RT(selectedRTIdx) ;
    std::vector<double> density(selectedRT.size()) ;

    if (!isempty[k])
    {
      if (bound[k])  // choose g_plus
      {
        for (arma::vec::iterator g_plusIt = selectedRT.begin();
          g_plusIt != selectedRT.end(); ++g_plusIt )
        {
          // a   v   zr   d  szr  sv  t0 st0 RT, precision;
          int l = std::distance(selectedRT.begin(), g_plusIt) ;
          gVec.at(8) = *g_plusIt ;
          if (checkDDM(gVec)) // Andrew's bad function
          {
            // 'cause the sampler may guess out-of-range pVec, catch it here and
            // set it to minLike density. This saves time for computing
            // expensive g_plus/g_minus function; should move it right
            // after the sampler to save more time
            out.at(selectedRTIdx[l]) = minLike ;
          }
          else
          {
            double tmp = std::abs(g_plus(gVec)) ;
            out.at(selectedRTIdx[l]) = std::max(tmp, minLike) ;
          }
        }
      }
      else
      {
        for (arma::vec::iterator g_minusIt = selectedRT.begin();
          g_minusIt != selectedRT.end(); ++g_minusIt )
        {
          // a   v   zr   d  szr  sv  t0 st0 RT, precision;
          int m = std::distance(selectedRT.begin(), g_minusIt) ;
          gVec.at(8) = *g_minusIt ;
          if (checkDDM(gVec))
          {
            out.at(selectedRTIdx[m]) = minLike ;
          }
          else
          {
            double tmp = std::abs(g_minus(gVec)) ;
            out.at(selectedRTIdx[m]) = std::max(tmp, minLike) ;
          }
        }
      }
    }

  }  // end of conditionNames loop (use Ctrl+p to track location)
  return out ;
}

//' @rdname ddmc
//' @export
// [[Rcpp::export]]
std::vector<double> ddmc_parallel(Rcpp::DataFrame x, Rcpp::NumericVector pVec,
  double precision=2.5, double minLike=1e-10)
{
  /* -------------------------------------------------------------------------
  First part - Create a boolean vector for choosing g_plus(1) or g_minus(0)
  ------------------------------------------------------------------------- */
  // Follow dnorm's sequence, I put quantile (named x) first (choice RTs) and
  // then the parameters (ie pVec); choice is achieved by ifelse "bound"
  // (see below). Others arguments eg fastdm's precision and GPU settings
  // come after (GPU setting has not added yet)
  Rcpp::NumericVector model = x.attr("model") ;         // For ddmc
  Rcpp::List matchMapOuter  = model.attr("match.map") ; // Handle two-layer list
  Rcpp::List matchMapInner  = matchMapOuter["M"] ;      // "M" as match.map's name
  Rcpp::List modelDimension = model.attr("dimnames") ;

  std::vector<std::string> stimulusName   = matchMapInner.names() ;
  std::vector<std::string> conditionNames = modelDimension[0] ; // row
  std::vector<std::string> parameterNames = modelDimension[1] ; // col; para
  std::vector<std::string> accumulatNames = modelDimension[2] ; // slice; r1, r2
  std::vector<int> bound (conditionNames.size());

  for (Rcpp::List::iterator it1 = matchMapInner.begin();
    it1 != matchMapInner.end(); ++it1)
  {
    int i = std::distance(matchMapInner.begin(), it1) ;
    std::string a = stimulusName[i] ;  // s1, s2
    std::string b = *it1 ;  // *it == matchMapInner[i]; r1, r2
    std::regex matchPattern (a + "(.*)" + b) ;
    for (std::vector<std::string>::iterator it2 = conditionNames.begin();
      it2 != conditionNames.end(); ++it2 )
    { // bound needs j
      int j = std::distance(conditionNames.begin(), it2) ;
      // Later accumulator test will overwrite previous, so test if only when
      // it is still not matched.
      if(bound[j] == 0)
      {
        bound[j] = (std::regex_match (*it2, matchPattern)) ? 1 : 0 ;
      }
    }
  }

  /* -------------------------------------------------------------------------
  Mawashigeri - min.like is set default in dmc.hpp
  ------------------------------------------------------------------------- */
  Rcpp::LogicalVector isempty = x.attr("cell.empty") ;
  Rcpp::List cellIdx = x.attr("cell.index") ; // s1.r1, s1.r2, s2.r1, s2.r2 etc.
  Rcpp::CharacterVector isemptyNames = isempty.names() ;
  arma::vec RT = x["RT"] ; // extract RTs from the input DataFrame (quantile)
  Rcpp::NumericMatrix accMat ; // forward declaration, avoid inside
  Rcpp::NumericVector accumulator ;     // loop.
  std::vector<double> out(RT.size()) ;  // ditto
  // prepare a container for pVec sending into g_minus/plus function, 'cause
  // my version of g_minus/plus require a redunctant element, RT
  std::vector<double> gVec(10) ;
  gVec.at(9) = precision ;
  arma::uvec selectedRTIdx ;  // 'cause arma's find function wants uvec type

  /* Let just use one vector, passing also RT and precision;
  * below is seq from dimnames(model)
  * 0  1  2  3  4   5   6  7    8   9
  * a  v  z  d  sz  sv  t0 st0  RT  precision (dmc uses names)
  */
  for (Rcpp::LogicalVector::iterator it3 = isempty.begin();
    it3 != isempty.end(); ++it3 )
  {
    // Why spend the effort to arrange the pVec for each accumulator and then
    // pick only the first accumulator?  Is it a bug? or
    // just to accommodate rtdists?
    int k = std::distance(isempty.begin(), it3) ; // n1order==true-> LBA
    accMat = getAccumulatorMatrix(pVec, conditionNames[k], model, true) ;
    accumulator = accMat(0, Rcpp::_) ;

    // fill in the gVec with the 8 parameter values firstly
    int gVecLen = gVec.size() - 2 ;
    for(int ii = 0; ii < gVecLen; ii++) { gVec[ii] = accumulator[ii] ; }
    // Because cellIdx is a List, I can use name indexing
    // Select RTs in conditionNames[k], eg s1.r1, condition. In order to use
    // indexing, cast Rcpp's NumericVector to Arma's vec
    // > attr(data.model,"cell.empty")
    //   s1.r1 s2.r1 s1.r2 s2.r2
    //   FALSE FALSE FALSE FALSE
    Rcpp::NumericVector cppIdxRT = cellIdx[conditionNames[k]] ; // T, T, F, etc.
    arma::vec armaIdxRT(cppIdxRT.begin(), cppIdxRT.size(), false) ;
    selectedRTIdx = find(armaIdxRT == 1) ;  // the location index equal TRUE
    arma::vec selectedRT = RT(selectedRTIdx) ;
    std::vector<double> density(selectedRT.size()) ;

    if (!isempty[k])
    {
      int nrt = selectedRT.size() ;
      if (bound[k])  // choose g_plus
      {
#pragma omp parallel for default(shared) firstprivate(gVec, selectedRT, selectedRTIdx)
        for (int iii=0; iii<nrt; iii++)
        {
          // a   v   z   d  sz  sv  t0 st0 RT, precision;
          gVec.at(8) = selectedRT[iii] ;
          if (checkDDM(gVec)) {out.at(selectedRTIdx[iii])=minLike ;}
          else
          {
            double tmp = std::abs(g_plus(gVec)) ;
            out.at(selectedRTIdx[iii]) = std::max(tmp, minLike) ;
          }
        }
      } else {
#pragma omp parallel for default(shared) firstprivate(gVec, selectedRT, selectedRTIdx)
        for (int jjj=0; jjj<nrt; jjj++)
        {
          gVec.at(8) = selectedRT[jjj] ;
          if (checkDDM(gVec)) {out.at(selectedRTIdx[jjj])=minLike ;}
          else
          {
            double tmp = std::abs(g_minus(gVec)) ;
            out.at(selectedRTIdx[jjj]) = std::max(tmp, minLike) ;
          }
        }
      }
    }

  }  // end of conditionNames loop (use Ctrl+p to track location)
  return out ;
}
