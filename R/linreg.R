# Documentation
#' A Reference Class to represent a linear regression package in R
#'
#' @description Listed functions are as methods in the object included:
#'
#' - print() prints out the coefficients and coefficient names, similar as done by the lm class
#'
#' - plot() plots the two plots using ggplot2
#'
#' 1. plot "Residuals vs Fitted"
#' 2. plot "Scale−Location"
#'
#' - resid() returns the vector of residuals e
#'
#' - pred() returns the predicted values y
#'
#' - coef() returns the coefficients as a named vector.
#'
#' - summary() returns a similar printout as printed for lm objects,
#'
#' - but only presents the coefficients with their standard error,
#'
#' - t-value and p-value as well as the estimate of σˆ and the degrees of freedom in the model.
#'
#'
#' @field X the values of X as matrix
#' @field y the values of y as matrix
#' @field regressions_coef the regression coefficients, saved as matrix
#' @field fitted_values the fitted values, saved as matrix
#' @field resi the residuals, saved as matrix
#' @field n the number of observations, saged as integer
#' @field p the number of parameters in the model, saved as integer
#' @field dof the degree of freedom, saved as integer
#' @field regressions_var the variance of the regression coefficent, saved as matrix
#' @field resi_var is the variance of residuals
#' @field t_value is the t-value, saved as matrix
#'
#' @field  m_formula is the formula
#' @field  m_data will return the name of the data
#'
#' @return Nothing
#' @export

linreg <- setRefClass("linreg",
                      fields = list(
                        X = "matrix",
                        y = "matrix",
                        regressions_coef = "matrix",
                        fitted_values = "matrix",
                        resi = "matrix",
                        n = "numeric",
                        p = "numeric",
                        dof = "numeric",
                        regressions_var = "matrix",
                        resi_var = "matrix",
                        t_value = "matrix",
                        m_formula ="formula",
                        m_data = "character"),

#-------------------------------------------------------------------------
                      methods = list
                      (
                        #This function will calculate necessary parameters
                        initialize = function(formula, data)
                        {
                          #Generate X and y
                          X <<- model.matrix(formula, data)
                          y <<- as.matrix(data[all.vars(formula)[1]])

                          #Regressions coeficients
                          regressions_coef<<- as.matrix(solve(t(X)%*%X) %*% t(X)%*%y)

                          fitted_values <<- X%*%regressions_coef

                          #Residuals
                          resi <<- y - fitted_values

                          #Degrees of fredom
                          n <<- length(X[,1])
                          p <<- length(X[1,])
                          dof <<- n - p

                          #Variance of the regression coeffcients
                          resi_var <<- (t(resi) %*% resi)/dof
                          regressions_var <<- as.numeric(resi_var) * solve(t(X) %*% X)

                          t_value <<- regressions_coef/as.double(sqrt(resi_var))

                          #Metadata
                          m_formula <<- formula
                          m_data <<- deparse(substitute(data))
                        },

                        #The print method
                        print = function() {
                          cat(paste("Call: \n"))
                          cat(paste("linreg(formula = ",format(m_formula), ", data = ", m_data, ")\n\n", sep = ""))
                          cat(paste("Coefficients:\n"))

                          coef <- structure(as.vector(regressions_coef), names= format(rownames(regressions_coef)))
                          my_print(coef)
                        },

                        #The plot method
                        plot = function()
                        {
                          library(ggplot2)

                          linkoping_theme <-
                            theme(
                              plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"),
                              panel.background = element_rect(fill="#3dd2dc"),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank(),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              axis.line = element_line(color= "#58585b", size=1),
                              axis.text.x = element_text(color="#1d1d1b", size="11"),
                              axis.text.y = element_text(color="#1d1d1b", size="11"),
                              axis.title.x = element_text(color="Black", size="12", face="bold.italic"),
                              axis.title.y = element_text(color="Black", size="12", face="bold.italic"),
                              axis.ticks.y = element_blank(),
                              axis.ticks.x = element_line(color = "#9D9ADC", size = 0.3),
                              plot.title = element_text(color="blue2", face="bold", size="14"),
                              legend.position="bottom", legend.title = element_blank(),
                              legend.text = element_text(color="Black", size="12")
                            )

                          #Plot 1
                          residuals_vs_fitted <- ggplot(data.frame(resi, fitted_values), aes(x=fitted_values, y=resi)) +
                            geom_point() +
                            stat_smooth(method='lm', colour="red", se=FALSE, span = 1) +
                            xlab(paste("Fitted Values\n", "linreg(", format(m_formula), ")", ""))+
                            ylab("Residuals") +
                            ggtitle("Residuals vs Fitted") +
                            linkoping_theme

                          #Print plot 1
                          my_print(residuals_vs_fitted)

                          #--------------------------------
                          #Prepare for plot 2
                          std_vs_fitted <- as.data.frame(cbind(sqrt(abs(resi-mean(resi))), fitted_values))
                          names(std_vs_fitted) = c("Standardized_residuals", "fitted_values")
                          y_plot <- std_vs_fitted[,1]

                          #Plot 2
                          plot2 <- ggplot(std_vs_fitted, aes(x = fitted_values, y = y_plot))+
                            geom_point() +
                            stat_smooth(method='lm', colour="red", se=FALSE, span = 1) +
                            xlab(paste("Fitted Values\n", "linreg(", format(m_formula), ")", ""))+
                            ylab(expression(sqrt("|Standardized residuals|"))) +
                            ggtitle("Scale Location") +
                            linkoping_theme

                          #Print plot2
                          my_print(plot2)
                        },
                        resid = function(){
                          return(resi)
                        },
                        pred = function()
                        {
                          return(fitted_values)
                        },
                        coef = function(){
                          coef <- structure(as.vector(regressions_coef), names= row.names(regressions_coef))
                          return(coef)
                        },

                        #Summary method
                        summary = function()
                        {
                          cat(paste("Call: \n\n"))
                          cat(paste("linreg(formula = ",format(m_formula), ", data = ", m_data, ")\n\n", sep = ""))
                          cat(paste("Coefficients:\n\n"))

                          regressions_var_diagonal <- diag(regressions_var)
                          table = data.frame(matrix(ncol = 5, nrow = 0))
                          for (i in 1:length(regressions_coef))
                          {
                            this_t_value = regressions_coef[i]/sqrt(regressions_var_diagonal[i])
                            this_p_value = 2*pt(abs(this_t_value), dof, lower.tail = FALSE)
                            row = data.frame(round(regressions_coef[i], 2), round(sqrt(regressions_var[i, i]), 2), round(this_t_value, 2), formatC(this_p_value, format = "e", digits = 2), write_star(this_p_value))
                            row.names(row)[1] = row.names(regressions_coef)[i]
                            table = rbind(table, row)
                          }
                          colnames(table) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
                          my_print(table)

                          cat(paste("\n"))
                          cat(paste("Residual standard error: ",round(sd(resi),3), " on ", dof, " degrees of freedom", sep = ""))
                        }

                      )
)

#'my_print function
#'Our own print function, because we can't call the default function inside RC object
#'@examples print("This is an example!")
#'@param x is the input.

my_print = function(x) {
    print(x)
}

#'write_star function
#'generate the start base on p_value, follow the rule of lm function
#'@param p_value the input
write_star = function(p_value) {
  if (p_value > 0.1)
    output <- ""
  else if (p_value > 0.05)
    output <- "."
  else if (p_value > 0.01)
    output <- "*"
  if (p_value > 0.001)
    output <- "**"
  else  output <- "***"
  return(output)
}

#How to use:
#1. Run this file (both linreg object and the my_print, write_star function)
#2. mod_object <- linreg(Petal.Length~Species, data = iris)
#3. mod_object$print()
#4. mod_object$plot()
#5. mod_object$summary()

#a <- linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data = iris)
