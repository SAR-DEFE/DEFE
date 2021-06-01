/*

Nombre:      ardl
Descripción: Do file para el documento de investigación "Pronósticos Uniecuacionales 
			 para los Ingresos Tributarios de Honduras". Publicado en la Revista de 
			 Administración Tributaria, N° 46, octubre 2020.
			 https://biblioteca.ciat.org/opac/book/5743
Autor:       Jose Carlo Bermúdez, jbermudez@sar.gob.hn

Departamento de Estudios Fiscales y Económicos
Direccion Nacional de Gestion Estratégica
Servicio de Administracion de Rentas, Honduras

*/

	clear all 
	clear matrix 
	global path " " 	// Ingrese directorio 
	cd "$path"
	use data_ardl.dta, clear

	tsmktim tiempo, start(2007m1) //Definicion de variable de tiempo y dummies estacionales (Definition of time series variable and seasonal dummies)
	tsset tiempo
	generate m=month(dofm(tiempo))
	tabulate m, generate(dum)
	gen trend=_n

	foreach var of varlist imae ener itrib imae_sa ener_sa itrib_sa {  //Rescalamiento de variables (Rescaling variables)
		gen ln_`var'=log(`var')
	}


	/*
	Gráfico de las variables
	*/
	tsline ln_imae ln_imae_sa if tiempo<=tm(2019m12),   ytitle("") xtitle("") ylabel(5(0.1)5.62) graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(blue red)
	tsline ln_itrib ln_itrib_sa if tiempo<=tm(2019m12), ytitle("") xtitle("") ylabel(7.7(0.3)9.53) graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(blue red)
	tsline ln_ener ln_ener_sa if tiempo<=tm(2019m12),   ytitle("") xtitle("") ylabel(12.8(0.1)13.45) graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(blue red)


	/*
	Rolling Correlations
	*/
	mvcorr ln_itrib_sa ln_imae_sa, win(12) gen(rho1) 
	mvcorr ln_itrib_sa ln_ener_sa, win(12) gen(rho2)

	bootstrap se_rho1=(r(sd)/sqrt(145)), reps(400): summ rho1 		//Error estándar mediante remuestreo (1) (Bootstrapped standar errors)
	bootstrap se_rho2=(r(sd)/sqrt(145)), reps(400): summ rho2

	gen rho1_upper=rho1+1.96*0.0301045				 				//Genero intervalos de confianza con errores estandar obtenidos con bootstrap (Confidence intervals with bootstrapped standard errors)
	gen rho1_lower=rho1-1.96*0.0201206
	gen rho2_upper=rho2+1.96*0.0304168
	gen rho2_lower=rho2-1.96*0.0240768

	graph twoway rarea rho1_upper rho1_lower tiempo if tiempo>=tm(2008m1), color(blue*.5) || line rho1 tiempo if tiempo>=tm(2008m1), lpattern(dash) lwidth(medthick) lcolor(blue*1.25)|| , ylabel(-1(0.5)1) ytitle("") xtitle("") graphregion(fcolor(white)) legend(row(2) ring(0) position(4) label(1 "± 95%") label(2 "Rho 1") order(2 1))
	graph twoway rarea rho2_upper rho2_lower tiempo if tiempo>=tm(2008m1), color(red*.5) || line rho2 tiempo if tiempo>=tm(2008m1), lpattern(dash) lwidth(medthick) lcolor(red*1.25)|| , ylabel(-1(0.5)1) ytitle("") xtitle("") graphregion(fcolor(white)) legend(row(2) ring(0) position(4) label(1 "± 95%") label(2 "Rho 2") order(2 1))


	/*
	Prueba de raíz unitaria empleando el algoritmo de Dolado et al (1990) // Unit root tests employing Dolado et al (1990) algorithm
	*/
	varsoc ln_imae_sa, exog(tiempo)                //Pruebas para el IMAE (IMAE unit root tests)
	dfuller ln_imae_sa, lags(3) trend reg
	dfuller ln_imae_sa, lags(3) reg
	dfuller ln_imae_sa, lags(3) noconstant reg

	varsoc ln_itrib_sa, exog(tiempo) 				//Pruebas para los Ingresos Tributarios (Tax revenue unit root tests)
	dfuller ln_itrib_sa, lags(2) trend reg
	dfuller ln_itrib_sa, lags(2) reg
	dfuller ln_itrib_sa, lags(2) noconstant reg

	varsoc ln_ener_sa, exog(tiempo)					//Pruebas para la Energia (Energy sales unit root tests)
	dfuller ln_ener_sa, lags(1) trend reg
	dfuller ln_ener_sa, lags(1) reg
	dfuller ln_ener_sa, lags(1) noconstant reg 


	/*
	Modelo ARDL mediante metodología de cointegración a la Pesaran, Shin & Smith (2001) \\ Estimation of the ARDL model by employing the cointegration approach of Pesaran et al (2001)
	*/
	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) aic dots  							//Estimacion de rezagos óptimos (Optimal lag estimation)
	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) ec1 lags(3 3 2) regstore(ardl_ec) 	//ARDL (p,q,r) en forma de MCE (Optimal lag model in its ECM representation)
	estat ectest 																									//Verifico Cointegración mediante prueba de Pesaran et al (2001) (checking for cointegration)
	estimates restore ardl_ec
	`e(cmdline)' vce(robust)
										
	predict ardl_1, residuals 						//Reservo los residuos del MCE (saving ECM residuals)
	estimates restore ardl_ec
	estat bgodfrey, lags(1) small					//Autocorrelación (checking for autocorrelation)
	estat hettest									//Homocedasticidad (checking for heterokedasticity)
	estat imtest, white

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl)			// Especificación Modelo ARDL Óptimo
	estimates restore ardl
	`e(cmdline)' vce(robust)


	/*
	Estimación de diferentes especificaciones de modelos univariados (Estimating different specifications of univariates models)
	*/
	forval p=1/12{
		qui arima ln_itrib, ar(`p') nolog //Modelo AR(p) óptimo
			estimates store ar_`p'
	}
	estimates stats ar_* 			 //Sugiere un AR(3)

	forval q=1/12{
		qui arima ln_itrib, ma(`q') nolog //Modelo MA(q) óptimo
			estimates store ma_`q'
	}
	estimates stats ma_* 			 //Sugiere un MA(8)

	forval p=1/4{
		forval q=1/4{
			qui arima ln_itrib, ar(`p') ma(`q') nolog //Modelo ARMA(p,q) óptimo
				estimates store arma_`p'_`q'
		}
	}
	estimates stats arma_*_*   			//Sugiere un ARMA(3,3)

	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib ln_imae if tin(, 2019m12), ar(1/`p') ma(1/`q') nolog //Modelo ARMAX(p,q) con el IMAE
				estimates store armaxi_`p'_`q'
		}
	}
	estimates stats armaxi_*_*			//Sugiere un ARMAX(2,1)

	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib ln_ener if tin(, 2019m12), ar(1/`p') ma(1/`q') nolog //Modelo ARMAX(p,q) con la Energía
				estimates store armaxe_`p'_`q'
		}
	}
	estimates stats armaxe_*_*			//Sugiere un ARMAX(1,1)

	set matsize 2000
	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib, arima(`p',1,0) sarima(`p',1,0,12) nolog //Modelo SAR(p,d,0)(P,D,Q)
				estimates store sar_`p'_1_0_`p'_1_0
			}
	}
	estimates stats sar_*_1_0_*_1_0		//Sugiere un SAR(3,0)(3,0)
		
	set matsize 2000
	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib, arima(`p',1,`q') sarima(`p',1,`q',12) nolog //Modelo SARIMA(p,d,q)(P,D,Q)
				estimates store sarima_`p'_1_`q'_`p'_1_`q'
			}
		}
	estimates stats sarima_*_1_*_*_1_*	//Sugiere un SARIMA(2,1)(2,1)


	/*
	Análisis de pronóstico dentro de muestra h=1 (One step ahead prediction capability of each model specification)
	*/
	qui arima ln_itrib, ar(3) vce(robust) 
	predict fitted_ar 

	qui arima ln_itrib, ma(8) vce(robust)
	predict fitted_ma

	qui arima ln_itrib, ar(3) ma(3) vce(robust)
	predict fitted_arma

	qui arima ln_itrib ln_imae, ar(2) ma(1) vce(robust)
	predict fitted_armaxi

	qui arima ln_itrib ln_ener, ar(1) ma(1) vce(robust)
	predict fitted_armaxe

	set matsize 2000
	qui arima ln_itrib, arima(3,1,0) sarima(3,1,0,12) vce(robust) nolog
	predict fitted_sar

	set matsize 2000
	qui arima ln_itrib, arima(2,1,1) sarima(2,1,1,12) vce(robust) nolog
	predict fitted_sarima

	forval p=1/3 {
		forval q=1/2 {
			quietly reg ln_itrib L(1/`p').ln_itrib L(1/`p').ln_ener L(1/`q').ln_imae dum1 dum4 dum6 dum8 dum10 dum12 trend, vce(robust)  //Modelos ARDL (ARDL models)
				predict ardl_fitted_`p'_`p'_`q'
		}
	}

	foreach var of varlist fitted_* ardl_fitted_*_*{
		fcstats ln_itrib `var'   								//Análisis de precisión pronósticos para h=1		
	}
			
	local fitted "ardl_fitted_3_3_2 fitted_ar fitted_ma fitted_arma fitted_armaxi fitted_armaxe ardl_fitted_1_1_1 ardl_fitted_1_1_2 ardl_fitted_2_2_1 ardl_fitted_2_2_2 ardl_fitted_3_3_1"
	tsline ln_itrib `fitted' if tiempo>=tm(2015m1), ylabel(8.2(0.5)9.8) ytitle("") xtitle("") graphregion(fcolor(white)) lwidth(medthick medthick thick thick thick thick thick thick thick thick thick thick) lcolor(blue red ltblue ltblue ltblue ltblue ltblue ltblue ltblue ltblue ltblue ltblue) lpattern(solid dash dot dot dot dot dot dot dot dot dot dot) legend(order(1 "Observado" 2 "ARDL(3,3,2)" 3 "Ajustados"))


	/*
	Análisis de pronóstico fuera de muestra h=3 (three periods ahead prediction capability of each model specification)
	*/
	qui arima ln_itrib if tiempo<=tm(2019m9), ar(3) vce(robust) 								// Especificación Modelo AR (AR model specification)
	estimates store ar_3h
	forecast create ar_model_3h, replace
	forecast estimates ar_3h
	forecast solve, suffix(_ar_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m9), ma(8) vce(robust)									// Especificación Modelo MA (MA model specification)
	estimates store ma_3h
	forecast create ma_model_3h, replace
	forecast estimates ma_3h
	forecast solve, suffix(_ma_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m9), ar(3) ma(3) vce(robust) 							// Especificación Modelo ARMA (ARMA model specification)
	estimates store arma_3h
	forecast create arma_model_3h, replace
	forecast estimates arma_3h
	forecast solve, suffix(_arma_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib ln_imae if tiempo<=tm(2019m9), ar(2) ma(1) vce(robust)					// Especificación Modelo ARMAX con IAME (ARMAX model specification with IMAE)
	estimates store armaxi_3h
	forecast create armaxi_model_3h, replace
	forecast estimates armaxi_3h
	forecast solve, suffix(_armaxi_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib ln_ener if tiempo<=tm(2019m9), ar(1) ma(1) vce(robust)      				// Especificación Modelo ARMAX con Energía (ARMAX model specification with energy sales)
	estimates store armaxe_3h
	forecast create armaxe_model_3h, replace
	forecast estimates armaxe_3h
	forecast solve, suffix(_armaxe_3h) begin(tm(2019m10)) log(off)

	set matsize 2000
	qui arima ln_itrib if tiempo<=tm(2019m9), arima(3,1,0) sarima(3,1,0,12) vce(robust) nolog 	// Especificación Modelo SAR (SAR model specification)
	estimates store sar_3h
	forecast create sar_model_3h, replace
	forecast estimates sar_3h, names(litrib)
	forecast solve, suffix(_sar_3h) begin(tm(2019m10)) log(off)

	set matsize 2000
	qui arima ln_itrib if tiempo<=tm(2019m9), arima(2,1,1) sarima(2,1,1,12) vce(robust) nolog 	// Especificación Modelo SARIMA (SARIMA model specification)
	estimates store sarima_3h
	forecast create sarima_model_3h, replace
	forecast estimates sarima_3h, names(lnitrib)
	forecast solve, suffix(_sarima_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl_1_1_1) // Especificación ARDL 1
	estimates restore ardl_1_1_1
	`e(cmdline)' vce(robust)
	estimates store ardl_1
	forecast create ardl_1_model_3h, replace
	forecast estimates ardl_1
	forecast solve, suffix(_ardl_1_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl_1_1_2) // Especificación ARDL 2
	estimates restore ardl_1_1_2
	`e(cmdline)' vce(robust)
	estimates store ardl_2
	forecast create ardl_2_model_3h, replace
	forecast estimates ardl_2
	forecast solve, suffix(_ardl_2_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl_2_2_1) // Especificación ARDL 3
	estimates restore ardl_2_2_1
	`e(cmdline)' vce(robust)
	estimates store ardl_3
	forecast create ardl_3_model_3h, replace
	forecast estimates ardl_3
	forecast solve, suffix(_ardl_3_3h) begin(tm(2019m10)) log(off)

	quietly ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl_2_2_2) // Especificación ARDL 4
	estimates restore ardl_2_2_2
	`e(cmdline)' vce(robust)
	estimates store ardl_4
	forecast create ardl_4_model_3h, replace
	forecast estimates ardl_4
	forecast solve, suffix(_ardl_4_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl_3_3_1) // Especificación ARDL 5
	estimates restore ardl_3_3_1
	`e(cmdline)' vce(robust)
	estimates store ardl_5
	forecast create ardl_5_model_3h, replace
	forecast estimates ardl_5
	forecast solve, suffix(_ardl_5_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl_3_3_2) // Especificación ARDL 5
	estimates restore ardl_3_3_2
	`e(cmdline)' vce(robust)
	estimates store ardl_6
	forecast create ardl_6_model_3h, replace
	forecast estimates ardl_6
	forecast solve, suffix(_ardl_6_3h) begin(tm(2019m10)) log(off)

	foreach var of varlist ln_itrib_*_* ln_itrib_*_*_* {
		fcstats ln_itrib `var'   													//Análisis de precisión pronósticos para h=3		
	}

	foreach var of varlist ln_itrib_*_* ln_itrib_*_*_* {
		drop `var'   																//Análisis de precisión pronósticos para h=3		
	}


	/*
	Análisis de pronóstico fuera de muestra h=6 (six periods ahead prediction capability of each model specification)
	*/
	qui arima ln_itrib if tiempo<=tm(2019m6), ar(3) vce(robust) 								// Especificación Modelo AR (AR model specification)
	estimates store ar_6h
	forecast create ar_model_6h, replace
	forecast estimates ar_6h
	forecast solve, suffix(_ar_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m6), ma(8) vce(robust)									// Especificación Modelo MA (MA model specification)
	estimates store ma_6h
	forecast create ma_model_6h, replace
	forecast estimates ma_6h
	forecast solve, suffix(_ma_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m6), ar(3) ma(3) vce(robust) 							// Especificación Modelo ARMA (ARMA model specification)
	estimates store arma_6h
	forecast create arma_model_6h, replace
	forecast estimates arma_6h
	forecast solve, suffix(_arma_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib ln_imae if tiempo<=tm(2019m6), ar(2) ma(1) vce(robust)					// Especificación Modelo ARMAX con IAME (ARMAX model specification with IMAE)
	estimates store armaxi_6h
	forecast create armaxi_model_6h, replace
	forecast estimates armaxi_6h
	forecast solve, suffix(_armaxi_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib ln_ener if tiempo<=tm(2019m6), ar(1) ma(1) vce(robust)      				// Especificación Modelo ARMAX con Energía (ARMAX model specification with energy sales)
	estimates store armaxe_6h
	forecast create armaxe_model_6h, replace
	forecast estimates armaxe_6h
	forecast solve, suffix(_armaxe_6h) begin(tm(2019m7)) log(off)

	set matsize 2000
	qui arima ln_itrib if tiempo<=tm(2019m6), arima(3,1,0) sarima(3,1,0,12) vce(robust) nolog 	// Especificación Modelo SAR (SAR model specification)
	estimates store sar_6h
	forecast create sar_model_6h, replace
	forecast estimates sar_6h, names(litrib1)
	forecast solve, suffix(_sar_6h) begin(tm(2019m7)) log(off)

	set matsize 2000
	qui arima ln_itrib if tiempo<=tm(2019m6), arima(2,1,1) sarima(2,1,1,12) vce(robust) nolog 	// Especificación Modelo SARIMA (SARIMA model specification)
	estimates store sarima_6h
	forecast create sarima_model_6h, replace
	forecast estimates sarima_6h, names(lnitrib1)
	forecast solve, suffix(_sarima_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl1_1_1_1) // Especificación ARDL 1
	estimates restore ardl1_1_1_1
	`e(cmdline)' vce(robust)
	estimates store ardl1_1
	forecast create ardl_1_model_6h, replace
	forecast estimates ardl1_1
	forecast solve, suffix(_ardl_1_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl1_1_1_2) // Especificación ARDL 2
	estimates restore ardl1_1_1_2
	`e(cmdline)' vce(robust)
	estimates store ardl1_2
	forecast create ardl_2_model_6h, replace
	forecast estimates ardl1_2
	forecast solve, suffix(_ardl_2_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl1_2_2_1) // Especificación ARDL 3
	estimates restore ardl1_2_2_1
	`e(cmdline)' vce(robust)
	estimates store ardl1_3
	forecast create ardl_3_model_6h, replace
	forecast estimates ardl1_3
	forecast solve, suffix(_ardl_3_6h) begin(tm(2019m7)) log(off)

	quietly ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl1_2_2_2) // Especificación ARDL 4
	estimates restore ardl1_2_2_2
	`e(cmdline)' vce(robust)
	estimates store ardl1_4
	forecast create ardl_4_model_6h, replace
	forecast estimates ardl1_4
	forecast solve, suffix(_ardl_4_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl1_3_3_1) // Especificación ARDL 5
	estimates restore ardl1_3_3_1
	`e(cmdline)' vce(robust)
	estimates store ardl1_5
	forecast create ardl_5_model_6h, replace
	forecast estimates ardl1_5
	forecast solve, suffix(_ardl_5_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl1_3_3_2) // Especificación ARDL 5
	estimates restore ardl1_3_3_2
	`e(cmdline)' vce(robust)
	estimates store ardl1_6
	forecast create ardl_6_model_6h, replace
	forecast estimates ardl1_6
	forecast solve, suffix(_ardl_6_6h) begin(tm(2019m7)) log(off)

	foreach var of varlist ln_itrib_*_* ln_itrib_*_*_* {
		fcstats ln_itrib `var'   													//Análisis de precisión pronósticos para h=6		
	}

	foreach var of varlist ln_itrib_*_* ln_itrib_*_*_* {
		drop `var'   																//Análisis de precisión pronósticos para h=6		
	}


	/*
	Análisis de pronóstico fuera de muestra h=12 (twelve periods ahead prediction capability of each model specification)
	*/
	qui arima ln_itrib if tiempo<=tm(2018m12), ar(3) vce(robust) 								// Especificación Modelo AR (AR model specification)
	estimates store ar_12h
	forecast create ar_model_12h, replace
	forecast estimates ar_12h
	forecast solve, suffix(_ar_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib if tiempo<=tm(2018m12), ma(8) vce(robust)									// Especificación Modelo MA (MA model specification)
	estimates store ma_12h
	forecast create ma_model_12h, replace
	forecast estimates ma_12h
	forecast solve, suffix(_ma_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib if tiempo<=tm(2018m12), ar(3) ma(3) vce(robust) 							// Especificación Modelo ARMA (ARMA model specification)
	estimates store arma_12h
	forecast create arma_model_12h, replace
	forecast estimates arma_12h
	forecast solve, suffix(_arma_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib ln_imae if tiempo<=tm(2018m12), ar(2) ma(1) vce(robust)					// Especificación Modelo ARMAX con IAME (ARMAX model specification with IMAE)
	estimates store armaxi_12h
	forecast create armaxi_model_12h, replace
	forecast estimates armaxi_12h
	forecast solve, suffix(_armaxi_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib ln_ener if tiempo<=tm(2018m12), ar(1) ma(1) vce(robust)      				// Especificación Modelo ARMAX con Energía (ARMAX model specification with energy sales)
	estimates store armaxe_12h
	forecast create armaxe_model_12h, replace
	forecast estimates armaxe_12h
	forecast solve, suffix(_armaxe_12h) begin(tm(2019m1)) log(off)

	set matsize 2000
	qui arima ln_itrib if tiempo<=tm(2018m12), arima(3,1,0) sarima(3,1,0,12) vce(robust) nolog 	// Especificación Modelo SAR (SAR model specification)
	estimates store sar_12h
	forecast create sar_model_12h, replace
	forecast estimates sar_12h, names(litrib2)
	forecast solve, suffix(_sar_12h) begin(tm(2019m1)) log(off)

	set matsize 2000
	qui arima ln_itrib if tiempo<=tm(2018m12), arima(2,1,1) sarima(2,1,1,12) vce(robust) nolog 	// Especificación Modelo SARIMA (SARIMA model specification)
	estimates store sarima_12h
	forecast create sarima_model_12h, replace
	forecast estimates sarima_12h, names(lnitrib2)
	forecast solve, suffix(_sarima_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl2_1_1_1) // Especificación ARDL 1
	estimates restore ardl2_1_1_1
	`e(cmdline)' vce(robust)
	estimates store ardl2_1
	forecast create ardl_1_model_12h, replace
	forecast estimates ardl2_1
	forecast solve, suffix(_ardl_1_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl2_1_1_2) // Especificación ARDL 2
	estimates restore ardl2_1_1_2
	`e(cmdline)' vce(robust)
	estimates store ardl2_2
	forecast create ardl_2_model_12h, replace
	forecast estimates ardl2_2
	forecast solve, suffix(_ardl_2_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl2_2_2_1) // Especificación ARDL 3
	estimates restore ardl2_2_2_1
	`e(cmdline)' vce(robust)
	estimates store ardl2_3
	forecast create ardl_3_model_12h, replace
	forecast estimates ardl2_3
	forecast solve, suffix(_ardl_3_12h) begin(tm(2019m1)) log(off)

	quietly ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl2_2_2_2) // Especificación ARDL 4
	estimates restore ardl2_2_2_2
	`e(cmdline)' vce(robust)
	estimates store ardl2_4
	forecast create ardl_4_model_12h, replace
	forecast estimates ardl2_4
	forecast solve, suffix(_ardl_4_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl2_3_3_1) // Especificación ARDL 5
	estimates restore ardl2_3_3_1
	`e(cmdline)' vce(robust)
	estimates store ardl2_5
	forecast create ardl_5_model_12h, replace
	forecast estimates ardl2_5
	forecast solve, suffix(_ardl_5_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl2_3_3_2) // Especificación ARDL 5
	estimates restore ardl2_3_3_2
	`e(cmdline)' vce(robust)
	estimates store ardl2_6
	forecast create ardl_6_model_12h, replace
	forecast estimates ardl2_6
	forecast solve, suffix(_ardl_6_12h) begin(tm(2019m1)) log(off)

	foreach var of varlist ln_itrib_*_* ln_itrib_*_*_* {
		fcstats ln_itrib `var'   													//Análisis de precisión pronósticos para h=12		
	}

	foreach var of varlist ln_itrib_*_* ln_itrib_*_*_* {
		drop `var'   																//Análisis de precisión pronósticos para h=12		
	}


	/*
	Estimación de multiplicadores de variables endógenas (Estimation of endogenous variables multipliers)
	*/
	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl1) // ARDL (1 1 1)
	estimates restore ardl1
	`e(cmdline)' vce(robust)
	gen mult_imae1=_b[ln_imae] in 1 												//Multiplicador del IMAE con ARDL(1 1 1)
	replace mult_imae1=L.mult_imae1*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
		forval c=3/12{
			replace mult_imae1=L.mult_imae1*_b[L1.ln_itrib] in `c'
	}
	gen mult_ener1=_b[ln_ener] in 1
	replace mult_ener1=L.mult_ener1*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				//Multiplicador de las ventas de energía con ARDL (1 1 1)
		forval c=3/12{
			replace mult_ener1=L.mult_ener1*_b[L1.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl2) // ARDL (1 1 2)
	estimates restore ardl2
	`e(cmdline)' vce(robust)
	gen mult_imae2=_b[ln_imae] in 1 												//Multiplicador del IMAE con ARDL (1 1 2)
	replace mult_imae2=L.mult_imae2*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae2=L.mult_imae2*_b[L1.ln_itrib]+_b[L2.ln_imae] in 3
		forval c=4/12{
			replace mult_imae2=L.mult_imae2*_b[L1.ln_itrib] in `c'
	}
	gen mult_ener2=_b[ln_ener] in 1
	replace mult_ener2=L.mult_ener2*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				//Multiplicador de las ventas de energía con ARDL (1 1 2)
		forval c=3/12{
			replace mult_ener2=L.mult_ener2*_b[L1.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl3) // ARDL (2 2 1)
	estimates restore ardl3
	`e(cmdline)' vce(robust)
	gen mult_imae3=_b[ln_imae] in 1 												//Multiplicador del IMAE con ARDL(2 2 1)
	replace mult_imae3=L.mult_imae3*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
		forval c=3/12{
			replace mult_imae3=L.mult_imae3*_b[L2.ln_itrib] in `c'
	}
	gen mult_ener3=_b[ln_ener] in 1
	replace mult_ener3=L.mult_ener3*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				//Multiplicador de las ventas de energía con ARDL (2 2 1)
	replace mult_ener3=L.mult_ener3*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
		forval c=4/12{
			replace mult_ener3=L.mult_ener3*_b[L2.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl4) // ARDL (2 2 2)
	estimates restore ardl4
	`e(cmdline)' vce(robust)
	gen mult_imae4=_b[ln_imae] in 1 												//Multiplicador del IMAE con ARDL(2 2 2)
	replace mult_imae4=L.mult_imae4*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae4=L.mult_imae4*_b[L2.ln_itrib]+_b[L2.ln_imae] in 3
		forval c=4/12{
			replace mult_imae4=L.mult_imae4*_b[L2.ln_itrib] in `c'
	}
	gen mult_ener4=_b[ln_ener] in 1
	replace mult_ener4=L.mult_ener4*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				//Multiplicador de las ventas de energía con ARDL (2 2 2)
	replace mult_ener4=L.mult_ener4*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
		forval c=4/12{
			replace mult_ener4=L.mult_ener4*_b[L2.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl5) // ARDL (3 3 1)
	estimates restore ardl5
	`e(cmdline)' vce(robust)
	gen mult_imae5=_b[ln_imae] in 1 												//Multiplicador del IMAE con ARDL (3 3 1)
	replace mult_imae5=L.mult_imae5*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae5=L.mult_imae5*_b[L2.ln_itrib] in 3
		forval c=4/12{
			replace mult_imae5=L.mult_imae5*_b[L3.ln_itrib] in `c'
	}
	gen mult_ener5=_b[ln_ener] in 1
	replace mult_ener5=L.mult_ener5*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				//Multiplicador de las ventas de energía con ARDL (3 3 1)
	replace mult_ener5=L.mult_ener5*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
	replace mult_ener5=L.mult_ener5*_b[L3.ln_itrib]+_b[L3.ln_ener] in 4
		forval c=5/12{
			replace mult_ener5=L.mult_ener5*_b[L3.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl6) // ARDL (3 3 2) Modelo Base
	estimates restore ardl6
	`e(cmdline)' vce(robust)
	gen mult_imae6=_b[ln_imae] in 1 												//Multiplicador del IMAE con ARDL (3 3 2)
	replace mult_imae6=L.mult_imae6*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae6=L.mult_imae6*_b[L2.ln_itrib]+_b[L2.ln_imae] in 3
		forval c=4/12{
			replace mult_imae6=L.mult_imae6*_b[L3.ln_itrib] in `c'
	}
	gen mult_ener6=_b[ln_ener] in 1
	replace mult_ener6=L.mult_ener6*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				//Multiplicador de las ventas de energía con ARDL (3 3 2)
	replace mult_ener6=L.mult_ener6*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
	replace mult_ener6=L.mult_ener6*_b[L3.ln_itrib]+_b[L3.ln_ener] in 4
		forval c=5/12{
			replace mult_ener6=L.mult_ener6*_b[L3.ln_itrib] in `c'
	}

	local mult_imae "mult_imae1 mult_imae2 mult_imae3 mult_imae4 mult_imae5 mult_imae6"  //Gráfica de multiplicadores del IMAE (plotting IMAE multipliers)
	egen mult_imae_mean=rowmean(`mult_imae')
	local y=0
	line mult_imae_mean mult_imae6 trend if trend<=12, ylabel(-.7(0.3).8) ytitle("") yline(`y') xlabel(1(1)12) xtitle("") graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(red blue) lpattern(longdash solid) legend(order(1 "Promedio" 2 "ARDL(3,3,2)"))

	local mult_ener "mult_ener1 mult_ener2 mult_ener3 mult_ener4 mult_ener5 mult_ener6"  //Gráfica de multiplicadores de la energía (plotting energy sales multipliers)
	egen mult_ener_mean=rowmean(`mult_ener')
	local y=0
	line mult_ener_mean mult_ener6 trend if trend<=12, ylabel(-.01(0.2).67) ytitle("") yline(`y') xlabel(1(1)12) xtitle("") graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(red blue) lpattern(longdash solid) legend(order(1 "Promedio" 2 "ARDL(3,3,2)"))
