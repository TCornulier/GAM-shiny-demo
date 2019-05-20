#### GAM - KNOT WEIGTH AS INPUT & CONSTANT NUMBER OF KNOTS

##### Function to procduce the SPLINE
# x takes the covariate values ("X")
# y takes the response values ("Y")
fit.spline<- function(x=x, y=y, knots.N= NULL, knots.vec= NULL, basis.coef){
	min.x<- min(x)
	max.x<- max(x)
	if(is.null(knots.vec)){
		knots.vec.std<- (1:knots.N) / (knots.N + 1)
	}
	else{
		knots.vec.std<- (knots.vec - min.x) / (max.x - min.x)
		if(!is.null(knots.N)){warning("Both knots.N and knots.vec have been specified. knots.N ignored.")}
	}
	if(length(basis.coef) != (length(knots.vec.std) + 2)){stop("'basis.coef' should be a vector of length equal to the number of knots + 2")}
	x.std<- (x - min.x) / (max.x - min.x)
	q<- length(knots.vec.std) + 2
	n<- length(x.std)
	X<- matrix(1, n, q) # fill with 1s. Component 1 = X[, 1] will be the intercept of the cubic spline
	X[,2]<- x.std # Component 2 = linear trend
	X[,3:q]<- 100 * outer(x.std, knots.vec.std, FUN= function(x, z){
		((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4 -((abs(x-z)-0.5)^4 -0.5*(abs(x-z)-0.5)^2 +7/240)/24 # computes the spline bases components of the cubic spline
	})
	weighted.bases<- t(t(X) * basis.coef) # multiply the cubic spline components by their respective basis coefficient
	spl.fit<- rowSums(weighted.bases) # fitted cubic spline is the sum of all components (spline bases)
	ESS= sum((x - spl.fit)^2) # Error Sum of Squares
    ### Output
    out.sp<-list(knots.vec.std=knots.vec.std, min.x=min.x, max.x=max.x, weighted.bases=weighted.bases, spl.fit=spl.fit, basis.coef=basis.coef,
				ESS=ESS)

    return(out.sp)
}


server<-function(input, output) {

    #### Input
	sp.data<-reactive({

    knots.N<-5

    basis.coef<-c(input$knot1, input$knot2, input$knot3, input$knot4, input$knot5, input$knot6, input$knot7)	#### Initial coeffs as inputs

    y<-sq$GSI[sq$Sex==2]

    x<-sq$MONTH[sq$Sex==2]

		fit.sp<-fit.spline(x=x, y=y, knots.N=knots.N, basis.coef=basis.coef)

		return(fit.sp)

		})


		output$caption<-renderText({

				ESS<-round(sp.data()$ESS, digits=2)				

				print(paste("Error Sum of Squares = ", ESS))

		})
  
    output$sp.plot<-renderPlot({

			### Data
			knots.N<-5

    	basis.coef<-sp.data()$basis.coef

    	y<-sq$GSI[sq$Sex==2]

    	x<-sq$MONTH[sq$Sex==2]


			### Plotting
      par(mfrow= c(1, 2))

			ord<-order(x)

			matplot(x[ord], sp.data()$weighted.bases[ord, ], type= "b", lty= 1, ylab= "Weighted basis functions", xlab= "Covariate (x)", main= "Cubic spline basis functions")

			abline(v= sp.data()$knots.vec.std * (sp.data()$max.x - sp.data()$min.x) + sp.data()$min.x, col= grey(0.8), lty= 1, lwd= 2)

	    		plot(y[ord] ~ x[ord], pch= "+", main= "Fitted Spline", xlab= "Covariate (x)", ylab= "Response (y)")

			abline(v=sp.data()$knots.vec.std * (sp.data()$max.x - sp.data()$min.x) + sp.data()$min.x, col= grey(0.8), lty= 1, lwd= 2)

			lines(x[ord], sp.data()$spl.fit[ord], lwd= 2)

    
        })
	

}
  
############ USER INTERFACE

ui<-fluidPage(theme=shinytheme("flatly"),

  titlePanel("Generalised Additive Model Example: 5 knots"),

	# Sidebar with controls to select the variable to plot against

    sidebarLayout(
    	sidebarPanel(

      		sliderInput("knot1",
                  "Weight of basis 1 (intercept):",
                  min = -10,
                  max = 10,
                  value = 1),

					sliderInput("knot2",
                  "Weight of basis 2 (linear trend):",
                  min = -10,
                  max = 10,
                  value = 1),
    
					sliderInput("knot3",
                  "Weight of basis 3 (spline basis for knot 1):",
                  min = -10,
                  max = 10,
                  value = 1),

				sliderInput("knot4",
                  "Weight of basis 4 (spline basis for knot 2):",
                  min = -10,
                  max = 10,
                  value = 1),

					sliderInput("knot5",
                  "Weight of basis 5 (spline basis for knot 3):",
                  min = -10,
                  max = 10,
                  value = 1),
					
					sliderInput("knot6",
                  "Weight of basis 6 (spline basis for knot 4):",
                  min = -10,
                  max = 10,
                  value = 1),

					sliderInput("knot7",
                  "Weight of basis 7 (spline basis for knot 5):",
                  min = -10,
                  max = 10,
                  value = 1)),

    mainPanel(
			h3(textOutput("caption")),
			plotOutput("sp.plot"))
  ),


	 h6("GAM, plots, and data: Thomas Cornulier"),
	 h6("Shiny App: Pablo Garcia Diaz")
)


shinyApp(ui=ui, server=server)
