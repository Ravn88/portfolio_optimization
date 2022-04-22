# [01] *** =======[ PACKAGES & LIBRARIES ]======= ***

  # Packages and libraries installation - STAGE I - Extraction, analysis and data visualisation
  
    # Data extraction
      install.packages("tidyr")
      install.packages("fpp")
      install.packages("dplyr")
    
    # Data analysis
      install.packages("fBasics")
    
    # Data visualisation
      install.packages("ggplot2")
      install.packages("corrplot")
      install.packages("plotly")
      install.packages("ggalt")
  
  
  # Packages and libraries installation - STAGE II - Modelling and advanced data analysis
  
    # LEVEL I - Price predictions - Stability and accuracy
      install.packages('rpart')
      install.packages('rattle')
      install.packages('rpart.plot')
      install.packages('RColorBrewer')
      install.packages('party')
      install.packages('MASS')
      install.packages('Hmisc')
      install.packages('splus2R')
    
    # LEVEL II - Portfolio optimization
      install.packages('DEoptim')
      install.packages('ROI')
      install.packages('tidyverse')
      install.packages('tidyquant')
      install.packages('PortfolioAnalytics')
      install.packages('PerformanceAnalytics')
      install.packages('zoo')
      install.packages('plotly')
      install.packages('timekit')
      install.packages('ggthemes')
      install.packages('timetk')
      install.packages('Rglpk')
  
  # Libraries loading
    library(tidyr)
    library(fpp)
    library(dplyr)
    
    library(fBasics)
    
    library(ggplot2)
    library(corrplot)
    library(plotly)
    library(ggalt)
    
    library(rpart)
    library(rattle)
    library(rpart.plot)
    library(RColorBrewer)
    library(party)
    library(MASS)
    library(Hmisc)
    library(splus2R)
    
    library(PortfolioAnalytics)
    library(PerformanceAnalytics)
    library(zoo)
    library(DEoptim)
    library(ROI)
    require(ROI.plugin.glpk)
    require(ROI.plugin.quadprog)
    library(tidyverse)
    library(tidyquant)
  

# [02] *** =======[ ETL - EXTRACTION, TRANSFORMATION, LOADING DATA ]======= ***
      
  # Data sets loading
  
    ETF_data <- read.csv("ETFs.csv", sep = ",", header = TRUE, na = c("N/A", ""))
    MUTUAL_data <- read.csv("Mutual Funds.csv", sep = ",", header = TRUE, na = c("N/A",""))
    
  
  # Data sets corrections
    
    ETF_data$price_earnings <- as.numeric(ETF_data$price_earnings)
    MUTUAL_data$price_earnings <- as.numeric(MUTUAL_data$price_earnings)
    
    ETF_data$fund_treynor_ratio_3years <- as.numeric(ETF_data$fund_treynor_ratio_3years)
    MUTUAL_data$fund_treynor_ratio_3years <- as.numeric(MUTUAL_data$fund_treynor_ratio_3years)
    
    MUTUAL_data$price_cashflow <- as.numeric(MUTUAL_data$price_cashflow)

    MUTUAL_data$price_sales <- as.numeric(MUTUAL_data$price_sales)
    
  
  # Removing columns from ETF data set which are not present in the MUTUAL FUNDS data set
    ETF_data$legal_type <- NULL
    
    
  # Removing columns from MUTUAL FUNDS data set which are not present in the ETF data set

    MUTUAL_data$morningstar_rating <- NULL
    MUTUAL_data$inception_date <- NULL
    MUTUAL_data$portfolio_cash <- NULL
    MUTUAL_data$portfolio_others <- NULL
    MUTUAL_data$portfolio_preferred <- NULL
    MUTUAL_data$portfolio_convertable <- NULL
    MUTUAL_data$median_market_cap <- NULL
    MUTUAL_data$bond_maturity <- NULL
    MUTUAL_data$bond_duration <- NULL
    MUTUAL_data$morningstar_return_rating <- NULL
    MUTUAL_data$category_return_2018 <- NULL
    MUTUAL_data$category_return_2017 <- NULL
    MUTUAL_data$category_return_2016 <- NULL
    MUTUAL_data$category_return_2015 <- NULL
    MUTUAL_data$category_return_2014 <- NULL
    MUTUAL_data$category_return_2013 <- NULL
    MUTUAL_data$category_return_2012 <- NULL
    MUTUAL_data$category_return_2011 <- NULL
    MUTUAL_data$category_return_2010 <- NULL
    MUTUAL_data$morningstar_risk_rating <- NULL
    MUTUAL_data$years_up <- NULL
    MUTUAL_data$years_down <- NULL
    
    
  # [OPTIONAL] - Removing columns from both ETF and MUTUAL data set which are common, but not necessary
    
    # From ETF
      ETF_data$fund_symbol <- NULL
      ETF_data$fund_name <- NULL
      ETF_data$fund_extended_name <- NULL
      ETF_data$category <- NULL
      ETF_data$fund_family <- NULL
      ETF_data$currency <- NULL
      ETF_data$rating_us_government <- NULL
      ETF_data$returns_ytd <- NULL
      ETF_data$alpha_3y <- NULL
      ETF_data$beta_3y <- NULL
      ETF_data$treynor_ratio_3y <- NULL
      ETF_data$sharpe_ratio_3y <- NULL
      ETF_data$standard_deviation_3y <- NULL
      ETF_data$mean_annual_return_3y <- NULL
      ETF_data$net_assets <- NULL
      ETF_data$ytd_return <- NULL
      ETF_data$fund_yield <- NULL
      ETF_data$inception_date <- NULL
      ETF_data$currency <- NULL
      ETF_data$net_annual_expense_ratio_fund <- NULL
      ETF_data$net_annual_expense_ratio_category <- NULL
      ETF_data$price_book <- NULL
      ETF_data$price_sales <- NULL
      ETF_data$price_cashflow <- NULL
      ETF_data$rating_us_government <- NULL
      ETF_data$fund_return_ytd <- NULL
      ETF_data$category_return_ytd <- NULL
      ETF_data$fund_return_1month <- NULL
      ETF_data$category_return_1month <- NULL
      ETF_data$fund_return_3months <- NULL
      ETF_data$category_return_3months <- NULL
      ETF_data$fund_return_1year <- NULL
      ETF_data$category_return_1year <- NULL
      ETF_data$fund_return_3years <- NULL
      ETF_data$category_return_3years <- NULL
      ETF_data$fund_return_5years <- NULL
      ETF_data$category_return_5years <- NULL
      ETF_data$fund_return_10years <- NULL
      ETF_data$category_return_10years <- NULL
      ETF_data$fund_alpha_3years <- NULL
      ETF_data$category_alpha_3years <- NULL
      ETF_data$fund_alpha_5years <- NULL
      ETF_data$category_alpha_5years <- NULL
      ETF_data$fund_alpha_10years <- NULL
      ETF_data$category_alpha_10years <- NULL
      ETF_data$fund_beta_3years <- NULL
      ETF_data$category_beta_3years <- NULL
      ETF_data$fund_beta_5years <- NULL
      ETF_data$category_beta_5years <- NULL
      ETF_data$fund_beta_10years <- NULL
      ETF_data$category_beta_10years <- NULL
      ETF_data$fund_mean_annual_return_3years <- NULL
      ETF_data$category_mean_annual_return_3years <- NULL
      ETF_data$fund_mean_annual_return_5years <- NULL
      ETF_data$category_mean_annual_return_5years <- NULL
      ETF_data$fund_mean_annual_return_10years <- NULL
      ETF_data$category_mean_annual_return_10years <- NULL
      ETF_data$fund_r_squared_3years <- NULL
      ETF_data$category_r_squared_3years <- NULL
      ETF_data$fund_r_squared_5years <- NULL
      ETF_data$category_r_squared_5years <- NULL
      ETF_data$fund_r_squared_10years <- NULL
      ETF_data$category_r_squared_10years <- NULL
      ETF_data$fund_standard_deviation_3years <- NULL
      ETF_data$category_standard_deviation_3years <- NULL
      ETF_data$fund_standard_deviation_5years <- NULL
      ETF_data$category_standard_deviation_5years <- NULL
      ETF_data$fund_standard_deviation_10years <- NULL
      ETF_data$category_standard_deviation_10years <- NULL
      ETF_data$fund_sharpe_ratio_3years <- NULL
      ETF_data$category_sharpe_ratio_3years <- NULL
      ETF_data$fund_sharpe_ratio_5years <- NULL
      ETF_data$category_sharpe_ratio_5years <- NULL
      ETF_data$fund_sharpe_ratio_10years <- NULL
      ETF_data$category_sharpe_ratio_10years <- NULL
      ETF_data$fund_treynor_ratio_3years <- NULL
      ETF_data$category_treynor_ratio_3years <- NULL
      ETF_data$fund_treynor_ratio_5years <- NULL
      ETF_data$category_treynor_ratio_5years <- NULL
      ETF_data$fund_treynor_ratio_10years <- NULL
      ETF_data$category_treynor_ratio_10years <- NULL
    
  
    # From MUTUAL FUNDS
      MUTUAL_data$fund_symbol <- NULL
      MUTUAL_data$fund_name <- NULL
      MUTUAL_data$fund_extended_name <- NULL
      MUTUAL_data$category <- NULL
      MUTUAL_data$fund_family <- NULL
      MUTUAL_data$currency <- NULL
      MUTUAL_data$rating_us_government <- NULL
      MUTUAL_data$returns_ytd <- NULL
      MUTUAL_data$alpha_3y <- NULL
      MUTUAL_data$beta_3y <- NULL
      MUTUAL_data$treynor_ratio_3y <- NULL
      MUTUAL_data$sharpe_ratio_3y <- NULL
      MUTUAL_data$standard_deviation_3y <- NULL
      MUTUAL_data$mean_annual_return_3y <- NULL
      MUTUAL_data$net_assets <- NULL
      MUTUAL_data$ytd_return <- NULL
      MUTUAL_data$fund_yield <- NULL
      MUTUAL_data$inception_date <- NULL
      MUTUAL_data$currency <- NULL
      MUTUAL_data$net_annual_expense_ratio_fund <- NULL
      MUTUAL_data$net_annual_expense_ratio_category <- NULL
      MUTUAL_data$price_book <- NULL
      MUTUAL_data$price_sales <- NULL
      MUTUAL_data$price_cashflow <- NULL
      MUTUAL_data$rating_us_government <- NULL
      MUTUAL_data$fund_return_ytd <- NULL
      MUTUAL_data$category_return_ytd <- NULL
      MUTUAL_data$fund_return_1month <- NULL
      MUTUAL_data$category_return_1month <- NULL
      MUTUAL_data$fund_return_3months <- NULL
      MUTUAL_data$category_return_3months <- NULL
      MUTUAL_data$fund_return_1year <- NULL
      MUTUAL_data$category_return_1year <- NULL
      MUTUAL_data$fund_return_3years <- NULL
      MUTUAL_data$category_return_3years <- NULL
      MUTUAL_data$fund_return_5years <- NULL
      MUTUAL_data$category_return_5years <- NULL
      MUTUAL_data$fund_return_10years <- NULL
      MUTUAL_data$category_return_10years <- NULL
      MUTUAL_data$fund_alpha_3years <- NULL
      MUTUAL_data$category_alpha_3years <- NULL
      MUTUAL_data$fund_alpha_5years <- NULL
      MUTUAL_data$category_alpha_5years <- NULL
      MUTUAL_data$fund_alpha_10years <- NULL
      MUTUAL_data$category_alpha_10years <- NULL
      MUTUAL_data$fund_beta_3years <- NULL
      MUTUAL_data$category_beta_3years <- NULL
      MUTUAL_data$fund_beta_5years <- NULL
      MUTUAL_data$category_beta_5years <- NULL
      MUTUAL_data$fund_beta_10years <- NULL
      MUTUAL_data$category_beta_10years <- NULL
      MUTUAL_data$fund_mean_annual_return_3years <- NULL
      MUTUAL_data$category_mean_annual_return_3years <- NULL
      MUTUAL_data$fund_mean_annual_return_5years <- NULL
      MUTUAL_data$category_mean_annual_return_5years <- NULL
      MUTUAL_data$fund_mean_annual_return_10years <- NULL
      MUTUAL_data$category_mean_annual_return_10years <- NULL
      MUTUAL_data$fund_r_squared_3years <- NULL
      MUTUAL_data$category_r_squared_3years <- NULL
      MUTUAL_data$fund_r_squared_5years <- NULL
      MUTUAL_data$category_r_squared_5years <- NULL
      MUTUAL_data$fund_r_squared_10years <- NULL
      MUTUAL_data$category_r_squared_10years <- NULL
      MUTUAL_data$fund_standard_deviation_3years <- NULL
      MUTUAL_data$category_standard_deviation_3years <- NULL
      MUTUAL_data$fund_standard_deviation_5years <- NULL
      MUTUAL_data$category_standard_deviation_5years <- NULL
      MUTUAL_data$fund_standard_deviation_10years <- NULL
      MUTUAL_data$category_standard_deviation_10years <- NULL
      MUTUAL_data$fund_sharpe_ratio_3years <- NULL
      MUTUAL_data$category_sharpe_ratio_3years <- NULL
      MUTUAL_data$fund_sharpe_ratio_5years <- NULL
      MUTUAL_data$category_sharpe_ratio_5years <- NULL
      MUTUAL_data$fund_sharpe_ratio_10years <- NULL
      MUTUAL_data$category_sharpe_ratio_10years <- NULL
      MUTUAL_data$fund_treynor_ratio_3years <- NULL
      MUTUAL_data$category_treynor_ratio_3years <- NULL
      MUTUAL_data$fund_treynor_ratio_5years <- NULL
      MUTUAL_data$category_treynor_ratio_5years <- NULL
      MUTUAL_data$fund_treynor_ratio_10years <- NULL
      MUTUAL_data$category_treynor_ratio_10years <- NULL
  
  
  # Removing data with "N/A" - there will be no row with missing data!
    ETF_data <- data.frame(na.omit(ETF_data))
    MUTUAL_data <- data.frame(na.omit(MUTUAL_data))
  
  # Display of all data in the console, for a quick check - before combining of both data sets
    print(ETF_data)
    View(ETF_data)
  
    print(MUTUAL_data)
    View(MUTUAL_data)

  # Combining of both data sets - ETF and Mutual Funds - in a new table
    BOTH_data <- rbind(ETF_data, MUTUAL_data)
    is.na(BOTH_data)
    
  # [OPTIONAL] - Removing the old data sets
    remove(ETF_data, MUTUAL_data)

  # Display of all data in the console, for a quick check - after combining of both data sets
    print(BOTH_data)
    View(BOTH_data)
  

# [03] *** =======[ VARIABLES MODIFYING AND DECLARATION ]======= ***

  # 1. Modification - converting asset classes weights and ratios into decimals (percents),
  # along assumption: 0 < weight_underclass <= 1 AND sum(weight_class) = 1
  
    # Weights and ratios of asset classes - bonds and equities
      BOTH_data$portfolio_bonds <- BOTH_data$portfolio_bonds / 100
      BOTH_data$portfolio_stocks <- BOTH_data$portfolio_stocks / 100
    
    # Converting return rates into decimals (percents)
      BOTH_data$fund_return_2010 <- BOTH_data$fund_return_2010 / 100
      BOTH_data$fund_return_2011 <- BOTH_data$fund_return_2011 / 100
      BOTH_data$fund_return_2012 <- BOTH_data$fund_return_2012 / 100
      BOTH_data$fund_return_2013 <- BOTH_data$fund_return_2013 / 100
      BOTH_data$fund_return_2014 <- BOTH_data$fund_return_2014 / 100
      BOTH_data$fund_return_2015 <- BOTH_data$fund_return_2015 / 100
      BOTH_data$fund_return_2016 <- BOTH_data$fund_return_2016 / 100
      BOTH_data$fund_return_2017 <- BOTH_data$fund_return_2017 / 100
      BOTH_data$fund_return_2018 <- BOTH_data$fund_return_2018 / 100
    
    # Converting sector weights within equities into decimals (percents)
      BOTH_data$basic_materials <- BOTH_data$basic_materials / 100
      BOTH_data$consumer_cyclical <- BOTH_data$consumer_cyclical / 100
      BOTH_data$financial_services <- BOTH_data$financial_services / 100
      BOTH_data$real_estate <- BOTH_data$real_estate / 100
      BOTH_data$consumer_defensive <- BOTH_data$consumer_defensive / 100
      BOTH_data$healthcare <- BOTH_data$healthcare / 100
      BOTH_data$utilities <- BOTH_data$utilities / 100
      BOTH_data$communication_services <- BOTH_data$communication_services / 100
      BOTH_data$energy <- BOTH_data$energy / 100
      BOTH_data$industrials <- BOTH_data$industrials / 100
      BOTH_data$technology <- BOTH_data$technology / 100
    
    # Converting rating weights within bonds into decimals (percents)
      BOTH_data$rating_aaa <- BOTH_data$rating_aaa / 100
      BOTH_data$rating_aa <- BOTH_data$rating_aa / 100
      BOTH_data$rating_a <- BOTH_data$rating_a / 100
      BOTH_data$rating_bbb <- BOTH_data$rating_bbb / 100
      BOTH_data$rating_bb <- BOTH_data$rating_bb / 100
      BOTH_data$rating_b <- BOTH_data$rating_b / 100
      BOTH_data$rating_below_b <- BOTH_data$rating_below_b / 100
      BOTH_data$rating_others <- BOTH_data$rating_others / 100
    
    
  # 2. Declaration - joining the weights of all assets classes in context of the weights of the whole portfolio,
  # along assumption: 0 < asset_weight <= 1 AND sum(all_assets_weight) = 1
  
    # Weights for each asset in context of the whole portfolio for equities, where:
    # asset_class(equities) * asset_underclass(equity_sector)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_basic_materials" = BOTH_data$basic_materials * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_consumer_cyclical" = BOTH_data$consumer_cyclical * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_financial_services" = BOTH_data$financial_services * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_real_estate" = BOTH_data$real_estate * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_consumer_defensive" = BOTH_data$consumer_defensive * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_healthcare" = BOTH_data$healthcare * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_utilities" = BOTH_data$utilities * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_communication_services" = BOTH_data$communication_services * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_energy" = BOTH_data$energy * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_industrials" = BOTH_data$industrials * BOTH_data$portfolio_stocks)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_STOCK_technology" = BOTH_data$technology * BOTH_data$portfolio_stocks)
    
    # Weights for each asset in context of the whole portfolio for bonds, where:
    # asset_class(bonds) * asset_underclass(bonds_rating)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_aaa" = BOTH_data$rating_aaa * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_aa" = BOTH_data$rating_aa * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_a" = BOTH_data$rating_a * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_bbb" = BOTH_data$rating_bbb * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_bb" = BOTH_data$rating_bb * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_b" = BOTH_data$rating_b * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_below_b" = BOTH_data$rating_below_b * BOTH_data$portfolio_bonds)
      BOTH_data <- cbind(BOTH_data, "WEIGHT_BOND_others" = BOTH_data$rating_others * BOTH_data$portfolio_bonds)
    
    # New weight not only for the whole class, but also for underclass of the assets which are not included in equities and bonds,
    # and which is a difference of the total for sum of weights and sum of equities and bonds weights
      BOTH_data <- cbind(BOTH_data, "WEIGHT_OTHER_assets" = (1 - (BOTH_data$portfolio_stocks + BOTH_data$portfolio_bonds)))
    

      
# [04] *** =======[ DATA ANALYSIS AND VISUALIZATION ]======= ***

  # Portfolio structure
        
    png(file = "01 - Equities proportions in the portfolios.jpg", width = 900, height = 900, units = "px")
    hist(
         BOTH_data$portfolio_stocks, 
      
         xlab="Equities ratio in portfolio", 
         ylab = "Portfolios quantity", 
         main = "Ratio of equities' share in portfolio", 
      
         col="forestgreen", 
         breaks = 10)
    dev.off()
    
    
    png(file = "02 - Bonds proportions in the portfolios.jpg", width = 900, height = 900, units = "px")
    hist(
         BOTH_data$portfolio_bonds, 
      
         xlab="Bonds ratio in portfolio", 
         ylab = "Portfolios quantity", 
         main = "Ratio of bonds' share in portfolio", 
      
         col="royalblue4", 
         breaks = 10)
    dev.off()
    
    
    png(file = "03 - Other assets proportions in the portfolios.jpg", width = 900, height = 900, units = "px")
    hist(
         BOTH_data$WEIGHT_OTHER_assets, 
      
         xlab="Other assets ratio in portfolio", 
         ylab = "Portfolios quantity", 
         main = "Ratio of other assets' share in portfolio", 
      
         col="khaki", 
         breaks = 10)
    dev.off()
    
    
    margines <- par(mar = c(4, 10.5, 1, 1) + 0.5)
    png(file = "04 - Sectors ratios.jpg", width = 900, height = 900, units = "px")
    boxplot(
            BOTH_data$basic_materials,
            BOTH_data$consumer_cyclical,
            BOTH_data$financial_services,
            BOTH_data$real_estate,
            BOTH_data$consumer_defensive,
            BOTH_data$healthcare,
            BOTH_data$utilities,
            BOTH_data$communication_services,
            BOTH_data$energy,
            BOTH_data$industrials,
            BOTH_data$technology,
            
            xlab = "Sector/equity ratio in portfolio", 
            main = "Distribution of sector/equity ratios in portfolio",
            names = c("Basic materials", "Consumer cyclical", "Financial services", "Real estate", "Consumer defensive", "Healthcare", "Utilities", "Communication services", "Energy", "Industrials", "Technology"),
            
            col = c("cadetblue3", "antiquewhite3", "aquamarine3", "chartreuse3", "azure3", "chocolate3", "bisque3", "mediumorchid3", "navajowhite3", "olivedrab3", "royalblue3"),
            las = 1,
            horizontal = TRUE)
    dev.off()
    
    
    margines <- par(mar = c(5, 5, 4, 2) + 0.1)
    png(file = "05 - Ratings ratios.jpg", width = 900, height = 900, units = "px")
    boxplot(
            BOTH_data$rating_aaa,
            BOTH_data$rating_aa,
            BOTH_data$rating_a,
            BOTH_data$rating_bbb,
            BOTH_data$rating_bb,
            BOTH_data$rating_b,
            BOTH_data$rating_below_b,
            BOTH_data$rating_others,
            
            xlab = "Rating/bond ratio in portfolio", 
            main = "Distribution of rating/bond ratios in portfolio",
            names = c("AAA", "AA", "A", "BBB", "BB", "B", "Below B", "Other"), 
            
            col = c("steelblue", "springgreen1", "springgreen2", "springgreen3", "springgreen4", "tomato1", "tomato2", "wheat3"),
            las=1, 
            horizontal = TRUE)
    dev.off()
    margines <- par(mar = c(4, 4, 4, 2) + 0.1)
  
  # Profits structure and distribution
    
    png(file = "06 - Price to earnings.jpg", width = 900, height = 900, units = "px")
    hist(
         BOTH_data$price_earnings, 
          
         xlab = "Price to earnings", 
         ylab = "Portfolios quantity", 
         main = "Distribution of price to earnings", 
          
         col = "lightgoldenrod3", 
         breaks = 10)
    dev.off()
    
    
    png(file = "07 - Return rates within years.jpg", width = 900, height = 900, units = "px")
    boxplot(
            BOTH_data$fund_return_2010,
            BOTH_data$fund_return_2011,
            BOTH_data$fund_return_2012,
            BOTH_data$fund_return_2013,
            BOTH_data$fund_return_2014,
            BOTH_data$fund_return_2015,
            BOTH_data$fund_return_2016,
            BOTH_data$fund_return_2017,
            BOTH_data$fund_return_2018,
            
            xlab = "Year", 
            ylab = "Return rate (%)", 
            names = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018"), 
            main = "Return rates in years 2010-2018",
            
            col = "turquoise")
    dev.off()
  
    
  # Correlation structure
    
    kolory_korelacji <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    korelacja_obligacje <- data.frame(
                                      BOTH_data$WEIGHT_BOND_aaa, 
                                      BOTH_data$WEIGHT_BOND_aa, 
                                      BOTH_data$WEIGHT_BOND_a, 
                                      BOTH_data$WEIGHT_BOND_bbb, 
                                      BOTH_data$WEIGHT_BOND_bb, 
                                      BOTH_data$WEIGHT_BOND_b, 
                                      BOTH_data$WEIGHT_BOND_below_b, 
                                      BOTH_data$WEIGHT_BOND_others)
    
    names(korelacja_obligacje)[1] <- "AAA"
    names(korelacja_obligacje)[2] <- "AA"
    names(korelacja_obligacje)[3] <- "A"
    names(korelacja_obligacje)[4] <- "BBB"
    names(korelacja_obligacje)[5] <- "BB"
    names(korelacja_obligacje)[6] <- "B"
    names(korelacja_obligacje)[7] <- "Below B"
    names(korelacja_obligacje)[8] <- "Other"
    
    
    png(file = "08 - Weight correlations of bonds' ratings within the portfolio.jpg", width = 900, height = 900, units = "px")
    corrplot(
             cor(korelacja_obligacje), 
             method ="color", 
             order = "original",
             col = kolory_korelacji(200), 
             type ="lower", 
             addCoef.col = "black", 
             tl.col ="black",
             tl.srt = 45, 
             title = "Weight correlations of bonds' ratings within the portfolio", 
             na.label = blank, 
             sig.level = 0.05, 
             insig = "blank", 
             mar = c(0, 0, 2, 0),
             diag = FALSE)
    dev.off()
    
    
    korelacja_akcje <- data.frame(
                                  BOTH_data$WEIGHT_STOCK_basic_materials, 
                                  BOTH_data$WEIGHT_STOCK_consumer_cyclical, 
                                  BOTH_data$WEIGHT_STOCK_financial_services, 
                                  BOTH_data$WEIGHT_STOCK_real_estate, 
                                  BOTH_data$WEIGHT_STOCK_consumer_defensive, 
                                  BOTH_data$WEIGHT_STOCK_healthcare, 
                                  BOTH_data$WEIGHT_STOCK_utilities, 
                                  BOTH_data$WEIGHT_STOCK_communication_services, 
                                  BOTH_data$WEIGHT_STOCK_energy, 
                                  BOTH_data$WEIGHT_STOCK_industrials, 
                                  BOTH_data$WEIGHT_STOCK_technology)
    
    names(korelacja_akcje)[1] <- "Basic materials"
    names(korelacja_akcje)[2] <- "Consumer cyclical"
    names(korelacja_akcje)[3] <- "Financial services"
    names(korelacja_akcje)[4] <- "Real estate"
    names(korelacja_akcje)[5] <- "Consumer defensive"
    names(korelacja_akcje)[6] <- "Healthcare"
    names(korelacja_akcje)[7] <- "Utilities"
    names(korelacja_akcje)[8] <- "Communication services"
    names(korelacja_akcje)[9] <- "Energy"
    names(korelacja_akcje)[10] <- "Industrials"
    names(korelacja_akcje)[11] <- "Technology"
    
    
    png(file = "09 - Weight correlations of equities' sectors within portfolio.jpg", width = 900, height = 900, units = "px")
    corrplot(
             cor(korelacja_akcje), 
             method ="color", 
             order = "original",
             col = kolory_korelacji(200), 
             type ="lower", 
             addCoef.col = "black", 
             tl.col = "black",
             tl.srt = 45, 
             title = "Weight correlations of equities' sectors within portfolio", 
             na.label = blank, 
             sig.level = 0.05, 
             insig = "blank", 
             mar = c(0, 0, 2, 0),
             diag = FALSE)
    dev.off()
    
    
    korelacja_komplet <- data.frame(
                                    korelacja_obligacje,
                                    korelacja_akcje,
                                    BOTH_data$WEIGHT_OTHER_assets)
    
    names(korelacja_komplet)[20] <- "Other assets"
    
    
    png(file = "10 - Weight correlations of all assets types within portfolio.jpg", width = 900, height = 900, units = "px")
    corrplot(
             cor(korelacja_komplet), 
             method = "color", 
             order = "original",
             col = kolory_korelacji(200), 
             type = "lower", 
             addCoef.col = "black", 
             tl.col = "black",
             tl.srt = 45, 
             title = "Weight correlations of all assets types within portfolio", 
             na.label = blank, 
             sig.level = 0.05, 
             insig = "blank", 
             mar = c(0, 0, 2, 0),
             diag = FALSE)
    dev.off()
    
    
    wartosci_kolowe <- c(sum(BOTH_data$investment == "Growth" & BOTH_data$size == "Small"), 
                         sum(BOTH_data$investment == "Growth" & BOTH_data$size == "Medium"), 
                         sum(BOTH_data$investment == "Growth" & BOTH_data$size == "Large"), 
                         sum(BOTH_data$investment == "Value" & BOTH_data$size == "Small"), 
                         sum(BOTH_data$investment == "Value" & BOTH_data$size == "Medium"), 
                         sum(BOTH_data$investment == "Value" & BOTH_data$size == "Large"), 
                         sum(BOTH_data$investment == "Blend" & BOTH_data$size == "Small"),
                         sum(BOTH_data$investment == "Blend" & BOTH_data$size == "Medium"),
                         sum(BOTH_data$investment == "Blend" & BOTH_data$size == "Large"))
    
    
    procenty_kolowe <- sprintf("%.2f%s", round(wartosci_kolowe*100 / sum(wartosci_kolowe), 2), "%")
    
    
    kolory_kolowe <- c("palevioletred1", "palevioletred2", "palevioletred3",
                       "palegreen1", "palegreen2", "palegreen3",
                       "lightskyblue1", "lightskyblue2", "lightskyblue3")
    
    
    etykiety_kolowe <- c("Growth\nSmall", 
                         "Growth\nMedium", 
                         "Growth\nLarge",
                         "Value\nSmall", 
                         "Value\nMedium", 
                         "Value\nLarge",
                         "Blend\nSmall", 
                         "Blend\nMedium", 
                         "Blend\nLarge")
    
    
    png(file = "11 - Portfolios by investment strategies and capitalization sizes.jpg", width = 900, height = 900, units = "px")
    
    pie(
        wartosci_kolowe, 
        labels = procenty_kolowe,
        main = "Portfolios by investment strategies/nand capitalization sizes", 
        col = kolory_kolowe,
        clockwise=TRUE,
        radius=0.7,
        border="white",
        cex = 1)
    
    legend("topright",
           etykiety_kolowe,
           cex = 1,
           bty = "n",
           fill = kolory_kolowe,
           y.intersp = 1.6)
    
    dev.off()


    
# [05] *** =======[ ADVANCED DATA ANALYSIS AND MODELLING ]======= ***
    
  # STAGE 01 - PRICE MODELLING AND PREDICTION

    # MODEL 01 - STABILITY - modelling of price stability in time due in context of evolution
    # with use of linear log
    
    
    # 01 - Installing and loading packages
    
      install.packages("caret", dependencies = TRUE)
      install.packager("mlbench")
      install.packages("kernlab")
      install.packages("randomForest")
      install.packages("LiblineaR")
      install.packages("ranger")
      install.packages("naivebayes")
      install.packages("kknn")
      install.packages("HiDimDA")
      install.packages("corrplot")
      
      library(caret)
      library(mlbench)
      library(kernlab)
      library(randomForest)
      library(LiblineaR)
      library(ranger)
      library(naivebayes)
      library(kknn)
      library(HiDimDA)
      library(corrplot)
    
    
    # 02 - Data Set
    
      portfolio_data_set <- data.frame(BOTH_data$fund_return_2010, 
                                       BOTH_data$fund_return_2011,
                                       BOTH_data$fund_return_2012,
                                       BOTH_data$fund_return_2013,
                                       BOTH_data$fund_return_2014,
                                       BOTH_data$fund_return_2015,
                                       BOTH_data$fund_return_2016,
                                       BOTH_data$fund_return_2017,
                                       BOTH_data$fund_return_2018)
      
      portfolio_data_set$fund_return_2018 <- as.factor(portfolio_data_set$fund_return_2018)
    
      
    # 03 - Setting Data Set for Training and Testing
    
      index_validation <- createDataPartition(y = portfolio_data_set[[1]], p=0.70, list=FALSE)
      
      TRAIN_70 <- portfolio_data_set[index_validation,]
      TEST_30 <-portfolio_data_set[-index_validation,]
    
    
    # 05 - Models training and testing
    
      control <- trainControl(method="cv", number=10)
      metric <- "RMSE"
      seed <- 7
    
      # Bagged CART (treebag)
        set.seed(seed)
        fit.treebag <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "treebag", metric = metric, trControl = control)
        
      # Bagged MARS (bagEarth)
        set.seed(seed)
        fit.bagEarth <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bagEarth", metric = metric, trControl = control)
        
      # Bagged MARS using gCV Pruning (bagEarthGCV)
        set.seed(seed)
        fit.bagEarthGCV <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bagEarthGCV", metric = metric, trControl = control)
        
      # Bayesian Generalized Linear Model (bayesglm)
        set.seed(seed)
        fit.bayesglm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bayesglm", metric = metric, trControl = control)
        
      # Boosted Generalized Additive Model (gamboost)
        set.seed(seed)
        fit.gamboost <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gamboost", metric = metric, trControl = control)
        
      # Boosted Generalized Linear Model (glmboost)
        set.seed(seed)
        fit.glmboost <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "glmboost", metric = metric, trControl = control)
        
      # Boosted Tree (blackboost)
        set.seed(seed)
        fit.blackboost <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "blackboost", metric = metric, trControl = control)
        
      # CART (rpart)
        set.seed(seed)
        fit.rpart <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rpart", metric = metric, trControl = control)
        
      # CART (rpart1SE)
        set.seed(seed)
        fit.rpart1SE <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rpart1SE", metric = metric, trControl = control)
        
      # CART (rpart2)
        set.seed(seed)
        fit.rpart2 <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rpart2", metric = metric, trControl = control)
        
      # Conditional Inference Random Forest (cforest)
        set.seed(seed)
        fit.cforest <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "cforest", metric = metric, trControl = control)
        
      # Conditional Inference Tree (ctree)
        set.seed(seed)
        fit.ctree <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "ctree", metric = metric, trControl = control)
          
      # Conditional Inference Tree (ctree2)
        set.seed(seed)
        fit.ctree2 <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "ctree2", metric = metric, trControl = control)
        
      # eXtreme Gradient Boosting (xgbDART)
        set.seed(seed)
        fit.xgbDART <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "xgbDART", metric = metric, trControl = control)
        
      # eXtreme Gradient Boosting (xgbTree)
        set.seed(seed)
        fit.xgbTree <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "xgbTree", metric = metric, trControl = control)
        
      # Generalized Linear Model (glm)
        set.seed(seed)
        fit.glm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "glm", metric = metric, trControl = control)
        
      # Generalized Linear Model with Stepwise Feature Selection (glmStepAIC)
        set.seed(seed)
        fit.glmStepAIC <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "glmStepAIC", metric = metric, trControl = control)
        
      # Linear Regression (lm)
        set.seed(seed)
        fit.lm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "lm", metric = metric, trControl = control)
        
      # Linear Regression with Stepwise Selection (lmStepAIC)
        set.seed(seed)
        fit.lmStepAIC <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "lmStepAIC", metric = metric, trControl = control)
        
      # Model Averaged Neural Network (avNNet)
        set.seed(seed)
        fit.avNNet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "avNNet", metric = metric, trControl = control)
        
      # Multivariate Adaptive Regression Spline (earth)
        set.seed(seed)
        fit.earth <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "earth", metric = metric, trControl = control)
        
      # Multivariate Adaptive Regression Splines (gcvEarth)
        set.seed(seed)
        fit.gcvEarth <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gcvEarth", metric = metric, trControl = control)
        
      # Negative Binomial Generalized Linear Model (glm.nb)
        set.seed(seed)
        fit.glm.nb <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "glm.nb", metric = metric, trControl = control)
        
      # Neural Network (nnet)
        set.seed(seed)
        fit.nnet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "nnet", metric = metric, trControl = control)
        
      # Neural Networks with Feature Extraction (pcaNNet)
        set.seed(seed)
        fit.pcaNNet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "pcaNNet", metric = metric, trControl = control)
        
      # Projection Pursuit Regression (ppr)
        set.seed(seed)
        fit.ppr <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "ppr", metric = metric, trControl = control)
        
      # Random Forest (ranger)
        set.seed(seed)
        fit.ranger <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "ranger", metric = metric, trControl = control)
        
      # Robust Linear Model (rlm)
        set.seed(seed)
        fit.rlm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rlm", metric = metric, trControl = control)
        
      # Stochastic Gradient Boosting (gbm)
        set.seed(seed)
        fit.gbm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gbm", metric = metric, trControl = control)
        
      # ZA DŁUGO
      # Tree Models from Genetic Algorithms (evtree)
        set.seed(seed)
        fit.evtree <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "evtree", metric = metric, trControl = control)
        
      # Bagged Logic Regression (logicBag)
        set.seed(seed)
        fit.logicBag <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "logicBag", metric = metric, trControl = control)
        
      # Bagged Model (bag)
        set.seed(seed)
        fit.bag <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bag", metric = metric, trControl = control)
        
      # Ensembles of Generalized Linear Models (randomGLM)
        set.seed(seed)
        fit.randomGLM <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "randomGLM", metric = metric, trControl = control)
        
      # Parallel Random Forest (parRF)
        set.seed(seed)
        fit.parRF <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "parRF", metric = metric, trControl = control)
        
      # Quantile Random Forest (qrf)
        set.seed(seed)
        fit.qrf <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "qrf", metric = metric, trControl = control)
        
      # Quantile Regression Neural Network (qrnn)
        set.seed(seed)
        fit.qrnn <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "qrnn", metric = metric, trControl = control)
        
      # Random Forest (Rborist)
        set.seed(seed)
        fit.Rborist <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "Rborist", metric = metric, trControl = control)
        
      # Random Forest (rf)
        set.seed(seed)
        fit.rf <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rf", metric = metric, trControl = control)
        
      # Random Forest by Randomization (extraTrees)
        set.seed(seed)
        fit.extraTrees <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "extraTrees", metric = metric, trControl = control)
        
      # Random Forest Rule-Based Model (rfRules)
        set.seed(seed)
        fit.rfRules <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rfRules", metric = metric, trControl = control)
        
      # Regularized Random Forest (RRF)
        set.seed(seed)
        fit.RRF <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "RRF", metric = metric, trControl = control)
        
      # Regularized Random Forest (RRFglobal)
        set.seed(seed)
        fit.RRFglobal <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "RRFglobal", metric = metric, trControl = control)
        
      # Bayesian Additive Regression Trees (bartMachine)
        set.seed(seed)
        fit.bartMachine <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bartMachine", metric = metric, trControl = control)
        
      # Bayesian Regularized Neural Networks (brnn)
        set.seed(seed)
        fit.brnn <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "brnn", metric = metric, trControl = control)
        
      # Bayesian Ridge Regression (bridge)
        set.seed(seed)
        fit.bridge <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bridge", metric = metric, trControl = control)
        
      # Bayesian Ridge Regression (Model Averaged) (blassoAveraged)
        set.seed(seed)
        fit.blassoAveraged <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "blassoAveraged", metric = metric, trControl = control)
        
      # Spike and Slab Regression (spikeslab)
        set.seed(seed)
        fit.spikeslab <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "spikeslab", metric = metric, trControl = control)
        
      # The Bayesian lasso (blasso)
        set.seed(seed)
        fit.blasso <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "blasso", metric = metric, trControl = control)
        
      # Logic Regression (logreg)
        set.seed(seed)
        fit.logreg <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "logreg", metric = metric, trControl = control)
        
      # Boosted Linear Model (BstLm)
        set.seed(seed)
        fit.BstLm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "BstLm", metric = metric, trControl = control)
        
      # Boosted Smoothing Spline (bstSm)
        set.seed(seed)
        fit.bstSm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bstSm", metric = metric, trControl = control)
        
      # Boosted Tree (bstTree)
        set.seed(seed)
        fit.bstTree <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bstTree", metric = metric, trControl = control)
        
      # Cubist (cubist)
        set.seed(seed)
        fit.cubist <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "cubist", metric = metric, trControl = control)
        
      # eXtreme Gradient Boosting (xgbLinear)
        set.seed(seed)
        fit.xgbLinear <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "xgbLinear", metric = metric, trControl = control)
        
      # Gradient Boosting Machines (gbm_h2o)
        set.seed(seed)
        fit.gbm_h2o <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gbm_h2o", metric = metric, trControl = control)
        
      # Tree-Based Ensembles (nodeHarvest)
        set.seed(seed)
        fit.nodeHarvest <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "nodeHarvest", metric = metric, trControl = control)
        
      # Independent Component Regression (icr)
        set.seed(seed)
        fit.icr <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "icr", metric = metric, trControl = control)
        
      # Partial Least Squares (kernelpls)
        set.seed(seed)
        fit.kernelpls <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "kernelpls", metric = metric, trControl = control)
        
      # Partial Least Squares (pls)
        set.seed(seed)
        fit.pls <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "pls", metric = metric, trControl = control)
        
      # Partial Least Squares (simpls)
        set.seed(seed)
        fit.simpls <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "simpls", metric = metric, trControl = control)
        
      # Partial Least Squares (widekernelpls)
        set.seed(seed)
        fit.widekernelpls <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "widekernelpls", metric = metric, trControl = control)
        
      # Principal Component Analysis (pcr)
        set.seed(seed)
        fit.pcr <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "pcr", metric = metric, trControl = control)
        
      # Sparse Partial Least Squares (spls)
        set.seed(seed)
        fit.spls <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "spls", metric = metric, trControl = control)
        
      # Supervised Principal Component Analysis (superpc)
        set.seed(seed)
        fit.superpc <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "superpc", metric = metric, trControl = control)
        
      # Linear Regression with Backwards Selection (leapBackward)
        set.seed(seed)
        fit.leapBackward <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "leapBackward", metric = metric, trControl = control)
        
      # Linear Regression with Forward Selection (leapForward)
        set.seed(seed)
        fit.leapForward <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "leapForward", metric = metric, trControl = control)
        
      # Linear Regression with Stepwise Selection (leapSeq)
        set.seed(seed)
        fit.leapSeq <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "leapSeq", metric = metric, trControl = control)
        
      # Ridge Regression with Variable Selection (foba)
        set.seed(seed)
        fit.foba <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "foba", metric = metric, trControl = control)
        
      # Gaussian Process (gaussprLinear)
        set.seed(seed)
        fit.gaussprLinear <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gaussprLinear", metric = metric, trControl = control)
        
      # Gaussian Process with Polynomial Kernel (gaussprPoly)
        set.seed(seed)
        fit.gaussprPoly <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gaussprPoly", metric = metric, trControl = control)
        
      # Gaussian Process with Radial Basis Function Kernel (gaussprRadial)
        set.seed(seed)
        fit.gaussprRadial <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gaussprRadial", metric = metric, trControl = control)
        
      # Generalized Additive Model using LOESS (gamLoess)
        set.seed(seed)
        fit.gamLoess <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gamLoess", metric = metric, trControl = control)
        
      # Generalized Additive Model using Splines (bam)
        set.seed(seed)
        fit.bam <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "bam", metric = metric, trControl = control)
        
      # Generalized Additive Model using Splines (gam)
        set.seed(seed)
        fit.gam <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gam", metric = metric, trControl = control)
        
      # Generalized Additive Model using Splines (gamSpline)
        set.seed(seed)
        fit.gamSpline <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "gamSpline", metric = metric, trControl = control)
        
      # glmnet (glmnet_h2o)
        set.seed(seed)
        fit.glmnet_h2o <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "glmnet_h2o", metric = metric, trControl = control)
        
      # glmnet (glmnet)
        set.seed(seed)
        fit.glmnet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "glmnet", metric = metric, trControl = control)
        
      # Multi-Step Adaptive MCP-Net (msaenet)
        set.seed(seed)
        fit.msaenet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "msaenet", metric = metric, trControl = control)
        
      # Elasticnet (enet)
        set.seed(seed)
        fit.enet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "enet", metric = metric, trControl = control)
        
      # Least Angle Regression (lars)
        set.seed(seed)
        fit.lars <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "lars", metric = metric, trControl = control)
        
      # Least Angle Regression (lars2)
        set.seed(seed)
        fit.lars2 <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "lars2", metric = metric, trControl = control)
        
      # Model Rules (M5Rules)
        set.seed(seed)
        fit.M5Rules <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "M5Rules", metric = metric, trControl = control)
        
      # Model Tree (M5)
        set.seed(seed)
        fit.M5 <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "M5", metric = metric, trControl = control)
        
      # Non-Convex Penalized Quantile Regression (rqnc)
        set.seed(seed)
        fit.rqnc <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rqnc", metric = metric, trControl = control)
        
      # Penalized Linear Regression (penalized)
        set.seed(seed)
        fit.penalized <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "penalized", metric = metric, trControl = control)
        
      # Quantile Regression with LASSO penalty (rqlasso)
        set.seed(seed)
        fit.rqlasso <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rqlasso", metric = metric, trControl = control)
        
      # Relaxed Lasso (relaxo)
        set.seed(seed)
        fit.relaxo <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "relaxo", metric = metric, trControl = control)
        
      # The lasso (lasso)
        set.seed(seed)
        fit.lasso <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "lasso", metric = metric, trControl = control)
        
      # L2 Regularized Support Vector Machine (dual) with Linear Kernel (svmLinear3)
        set.seed(seed)
        fit.svmLinear3 <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmLinear3", metric = metric, trControl = control)
        
      # Polynomial Kernel Regularized Least Squares (krlsPoly)
        set.seed(seed)
        fit.krlsPoly <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "krlsPoly", metric = metric, trControl = control)
        
      # Radial Basis Function Kernel Regularized Least Squares (krlsRadial)
        set.seed(seed)
        fit.krlsRadial <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "krlsRadial", metric = metric, trControl = control)
        
      # Relevance Vector Machines with Linear Kernel (rvmLinear)
        set.seed(seed)
        fit.rvmLinear <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rvmLinear", metric = metric, trControl = control)
        
      # Relevance Vector Machines with Polynomial Kernel (rvmPoly)
        set.seed(seed)
        fit.rvmPoly <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rvmPoly", metric = metric, trControl = control)
        
      # Relevance Vector Machines with Radial Basis Function Kernel (rvmRadial)
        set.seed(seed)
        fit.rvmRadial <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rvmRadial", metric = metric, trControl = control)
        
      # Support Vector Machines with Boundrange String Kernel (svmBoundrangeString)
        set.seed(seed)
        fit.svmBoundrangeString <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmBoundrangeString", metric = metric, trControl = control)
        
      # Support Vector Machines with Exponential String Kernel (svmExpoString)
        set.seed(seed)
        fit.svmExpoString <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmExpoString", metric = metric, trControl = control)
        
      # Support Vector Machines with Linear Kernel (svmLinear)
        set.seed(seed)
        fit.svmLinear <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmLinear", metric = metric, trControl = control)
        
      # Support Vector Machines with Linear Kernel (svmLinear2)
        set.seed(seed)
        fit.svmLinear2 <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmLinear2", metric = metric, trControl = control)
        
      # Support Vector Machines with Polynomial Kernel (svmPoly)
        set.seed(seed)
        fit.svmPoly <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmPoly", metric = metric, trControl = control)
        
      # Support Vector Machines with Radial Basis Function Kernel (svmRadial)
        set.seed(seed)
        fit.svmRadial <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmRadial", metric = metric, trControl = control)
        
      # Support Vector Machines with Radial Basis Function Kernel (svmRadialCost)
        set.seed(seed)
        fit.svmRadialCost <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmRadialCost", metric = metric, trControl = control)
        
      # Support Vector Machines with Radial Basis Function Kernel (svmRadialSigma)
        set.seed(seed)
        fit.svmRadialSigma <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmRadialSigma", metric = metric, trControl = control)
        
      # Support Vector Machines with Spectrum String Kernel (svmSpectrumString)
        set.seed(seed)
        fit.svmSpectrumString <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "svmSpectrumString", metric = metric, trControl = control)
        
      # Multi-Layer Perceptron (mlpWeightDecay)
        set.seed(seed)
        fit.mlpWeightDecay <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlpWeightDecay", metric = metric, trControl = control)
        
      # Multi-Layer Perceptron, multiple layers (mlpWeightDecayML)
        set.seed(seed)
        fit.mlpWeightDecayML <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlpWeightDecayML", metric = metric, trControl = control)
        
      # Multilayer Perceptron Network by Stochastic Gradient Descent (mlpSGD)
        set.seed(seed)
        fit.mlpSGD <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlpSGD", metric = metric, trControl = control)
        
      # Multilayer Perceptron Network with Weight Decay (mlpKerasDecay)
        set.seed(seed)
        fit.mlpKerasDecay <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlpKerasDecay", metric = metric, trControl = control)
        
      # Radial Basis Function Network (rbf)
        set.seed(seed)
        fit.rbf <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rbf", metric = metric, trControl = control)
        
      # Radial Basis Function Network (rbfDDA)
        set.seed(seed)
        fit.rbfDDA <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "rbfDDA", metric = metric, trControl = control)
        
      # Ridge Regression (ridge)
        set.seed(seed)
        fit.ridge <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "ridge", metric = metric, trControl = control)
        
      # Non-Negative Least Squares (nnls)
        set.seed(seed)
        fit.nnls <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "nnls", metric = metric, trControl = control)
        
      # Extreme Learning Machine (elm)
        set.seed(seed)
        fit.elm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "elm", metric = metric, trControl = control)
        
      # Monotone Multi-Layer Perceptron Neural Network (monmlp)
        set.seed(seed)
        fit.monmlp <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "monmlp", metric = metric, trControl = control)
        
      # Multi-Layer Perceptron (mlp)
        set.seed(seed)
        fit.mlp <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlp", metric = metric, trControl = control)
        
      # Multi-Layer Perceptron, with multiple layers (mlpML)
        set.seed(seed)
        fit.mlpML <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlpML", metric = metric, trControl = control)
        
      # Multilayer Perceptron Network with Dropout (mlpKerasDropout)
        set.seed(seed)
        fit.mlpKerasDropout <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mlpKerasDropout", metric = metric, trControl = control)
        
      # Neural Network (mxnet)
        set.seed(seed)
        fit.mxnet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mxnet", metric = metric, trControl = control)
        
      # Neural Network (mxnetAdam)
        set.seed(seed)
        fit.mxnetAdam <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "mxnetAdam", metric = metric, trControl = control)
        
      # Neural Network (neuralnet)
        set.seed(seed)
        fit.neuralnet <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "neuralnet", metric = metric, trControl = control)
        
      # Stacked AutoEncoder Deep Neural Network (dnn)
        set.seed(seed)
        fit.dnn <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "dnn", metric = metric, trControl = control)
        
      # Partial Least Squares Generalized Linear Models (plsRglm)
        set.seed(seed)
        fit.plsRglm <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "plsRglm", metric = metric, trControl = control)
        
      # k-Nearest Neighbors (kknn)
        set.seed(seed)
        fit.kknn <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "kknn", metric = metric, trControl = control)
        
      # k-Nearest Neighbors (knn)
        set.seed(seed)
        fit.knn <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "knn", metric = metric, trControl = control)
        
      # Adaptive-Network-Based Fuzzy Inference System (ANFIS)
        set.seed(seed)
        fit.ANFIS <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "ANFIS", metric = metric, trControl = control)
        
      # Dynamic Evolving Neural-Fuzzy Inference System (DENFIS)
        set.seed(seed)
        fit.DENFIS <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "DENFIS", metric = metric, trControl = control)
        
      # Fuzzy Inference Rules by Descent Method (FIR.DM)
        set.seed(seed)
        fit.FIR.DM <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "FIR.DM", metric = metric, trControl = control)
        
      # Fuzzy Rules via MOGUL (GFS.FR.MOGUL)
        set.seed(seed)
        fit.GFS.FR.MOGUL <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "GFS.FR.MOGUL", metric = metric, trControl = control)
        
      # Fuzzy Rules via Thrift (GFS.THRIFT)
        set.seed(seed)
        fit.GFS.THRIFT <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "GFS.THRIFT", metric = metric, trControl = control)
        
      # Genetic Lateral Tuning and Rule Selection of Linguistic Fuzzy Systems (GFS.LT.RS)
        set.seed(seed)
        fit.GFS.LT.RS <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "GFS.LT.RS", metric = metric, trControl = control)
        
      # Hybrid Neural Fuzzy Inference System (HYFIS)
        set.seed(seed)
        fit.HYFIS <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "HYFIS", metric = metric, trControl = control)
        
      # Simplified TSK Fuzzy Rules (FS.HGD)
        set.seed(seed)
        fit.FS.HGD <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "FS.HGD", metric = metric, trControl = control)
        
      # Subtractive Clustering and Fuzzy c-Means Rules (SBC)
        set.seed(seed)
        fit.SBC <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "SBC", metric = metric, trControl = control)
        
      # Wang and Mendel Fuzzy Rules (WM)
        set.seed(seed)
        fit.WM <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "WM", metric = metric, trControl = control)
        
      # Self-Organizing Maps (xyf)
        set.seed(seed)
        fit.xyf <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "xyf", metric = metric, trControl = control)
        
      # Non-Informative Model (null)
        set.seed(seed)
        fit.null <- train(BOTH_data.fund_return_2018~., data = TRAIN_70, method = "null", metric = metric, trControl = control)
    
    
    
    # 06 - Accuracy summary
    
      results <- resamples(list(lm = fit.lm, bridge = fit.bridge, treebag = fit.treebag))
      
      summary(results)
      
      # Accuracy comparison
        dotplot(results)
      
      # Best Model summary
        fit.best = fit.svmLinear2
        print(fit.best)
    
    
    # MODEL 02 - PREDYKCJA - modelowanie zwrot?w za 2019 rok na bazie danych z lat 2010-2018
    
      # Predictions and skill estimation of the best model on the testing dataset
      predictions <- predict(fit.best, TEST_30)
      confusionMatrix(predictions, TEST_30$fund_return_2018)
      
      
      
      
      # A - INITIAL RETURN RATE PREDICTION FOR EACH YEAR - in order to check if predicted return rates are stable in time
      # and due to portfolio optimization
      
      
        predyktory_2019 <- lm(BOTH_data$fund_return_2018 ~ 
                                
                                BOTH_data$fund_return_2017 + 
                                BOTH_data$fund_return_2016 + 
                                BOTH_data$fund_return_2015 + 
                                BOTH_data$fund_return_2014 + 
                                BOTH_data$fund_return_2013 + 
                                BOTH_data$fund_return_2012 + 
                                BOTH_data$fund_return_2011 + 
                                BOTH_data$fund_return_2010)
        
        BOTH_data <- cbind(BOTH_data, "predicted_return_2019" = predict(predyktory_2019))
      
      
      # B - JOINING NEW WEIGHTS - for each asset in the context of the whole portfolio
      # Returns for each asset class by its weight in the context of the whole portfolio
      
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_OTHER_ASSETS" = BOTH_data$WEIGHT_OTHER_assets * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_basic_materials" = BOTH_data$WEIGHT_STOCK_basic_materials * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_consumer_cyclical" = BOTH_data$WEIGHT_STOCK_consumer_cyclical * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_financial_services" = BOTH_data$WEIGHT_STOCK_financial_services * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_real_estate" = BOTH_data$WEIGHT_STOCK_real_estate * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_consumer_defensive" = BOTH_data$WEIGHT_STOCK_consumer_defensive * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_healthcare" = BOTH_data$WEIGHT_STOCK_healthcare * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_utilities" = BOTH_data$WEIGHT_STOCK_utilities * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_communication_services" = BOTH_data$WEIGHT_STOCK_communication_services * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_energy" = BOTH_data$WEIGHT_STOCK_energy * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_industrials" = BOTH_data$WEIGHT_STOCK_industrials * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_STOCK_technology" = BOTH_data$WEIGHT_STOCK_technology * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_aaa" = BOTH_data$WEIGHT_BOND_aaa * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_aa" = BOTH_data$WEIGHT_BOND_aa * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_a" = BOTH_data$WEIGHT_BOND_a * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_bbb" = BOTH_data$WEIGHT_BOND_bbb * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_bb" = BOTH_data$WEIGHT_BOND_bb * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_b" = BOTH_data$WEIGHT_BOND_b * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_below_b" = BOTH_data$WEIGHT_BOND_below_b * BOTH_data$predicted_return_2019)
        BOTH_data <- cbind(BOTH_data, "PREDICT_RETURN_BOND_others" = BOTH_data$WEIGHT_BOND_others * BOTH_data$predicted_return_2019)  
      
      
      
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  # STAGE 02 - PORTFOLIO OPTIMIZATION

    # MODEL 03 - OPTIMIZATION - modelling of optimal proportions of assets within portfolio due to return and risk
        
      # C - BUCKETS SEGMENTATION - creation of separate data sets based on investment strategy
        
        # BLEND investment strategy and portfolio sizes
          SET_BLEND_SMALL <- subset(BOTH_data, BOTH_data$investment == "Blend" & BOTH_data$size == "Small")
          SET_BLEND_MEDIUM <- subset(BOTH_data, BOTH_data$investment == "Blend" & BOTH_data$size == "Medium")
          SET_BLEND_LARGE <- subset(BOTH_data, BOTH_data$investment == "Blend" & BOTH_data$size == "Large")
        
        # VALUE investment strategy and portfolio sizes
          SET_VALUE_SMALL <- subset(BOTH_data, BOTH_data$investment == "Value" & BOTH_data$size == "Small")
          SET_VALUE_MEDIUM <- subset(BOTH_data, BOTH_data$investment == "Value" & BOTH_data$size == "Medium")
          SET_VALUE_LARGE <- subset(BOTH_data, BOTH_data$investment == "Value" & BOTH_data$size == "Large")
        
        # GROWTH investment strategy and portfolio sizes
          SET_GROWTH_SMALL <- subset(BOTH_data, BOTH_data$investment == "Growth" & BOTH_data$size == "Small")
          SET_GROWTH_MEDIUM <- subset(BOTH_data, BOTH_data$investment == "Growth" & BOTH_data$size == "Medium")
          SET_GROWTH_LARGE <- subset(BOTH_data, BOTH_data$investment == "Growth" & BOTH_data$size == "Large")  
        
        
      # STEP 0 - Transform of data set into time series
        
        TS_SET_BLEND_SMALL <- ts(SET_BLEND_SMALL, frequency = 1, start = c(1970,1))
        TS_SET_BLEND_MEDIUM <- ts(SET_BLEND_MEDIUM, frequency = 1, start = c(1970,1))
        TS_SET_BLEND_LARGE <- ts(SET_BLEND_LARGE, frequency = 1, start = c(1970,1))
        TS_SET_GROWTH_SMALL <- ts(SET_GROWTH_SMALL, frequency = 1, start = c(1970,1))
        TS_SET_GROWTH_MEDIUM <- ts(SET_GROWTH_MEDIUM, frequency = 1, start = c(1970,1))
        TS_SET_GROWTH_LARGE <- ts(SET_GROWTH_LARGE, frequency = 1, start = c(1970,1))
        TS_SET_VALUE_SMALL <- ts(SET_VALUE_SMALL, frequency = 1, start = c(1970,1))
        TS_SET_VALUE_MEDIUM <- ts(SET_VALUE_MEDIUM, frequency = 1, start = c(1970,1))
        TS_SET_VALUE_LARGE <- ts(SET_VALUE_LARGE, frequency = 1, start = c(1970,1))
      
        
      # STEP 1- Initial parametrization based on predicted returns
        
        BLEND_returns_SMALL <- TS_SET_BLEND_SMALL[, 55:74]
        BLEND_returns_MEDIUM <- TS_SET_BLEND_MEDIUM[, 55:74]
        BLEND_returns_LARGE <- TS_SET_BLEND_LARGE[, 55:74]
        GROWTH_returns_SMALL <- TS_SET_GROWTH_SMALL[, 55:74]
        GROWTH_returns_MEDIUM <- TS_SET_GROWTH_MEDIUM[, 55:74]
        GROWTH_returns_LARGE <- TS_SET_GROWTH_LARGE[, 55:74]
        VALUE_returns_SMALL <- TS_SET_VALUE_SMALL[, 55:74]
        VALUE_returns_MEDIUM <- TS_SET_VALUE_MEDIUM[, 55:74]
        VALUE_returns_LARGE <- TS_SET_VALUE_LARGE[, 55:74]
      
        lista_aktywow <- c("OTHER ASSETS",
                           
                           "EQUITIES - Basic materials", 
                           "EQUITIES - Consumer cyclical", 
                           "EQUITIES - Financial services", 
                           "EQUITIES - Real estate", 
                           "EQUITIES - Consumer defensive", 
                           "EQUITIES - Healthcare", 
                           "EQUITIES - Utilities", 
                           "EQUITIES - Communication services", 
                           "EQUITIES - Energy", 
                           "EQUITIES - Industrials", 
                           "EQUITIES - Technology", 
                           
                           "BONDS - AAA", 
                           "BONDS - AA", 
                           "BONDS - A", 
                           "BONDS - BBB", 
                           "BONDS - BB", 
                           "BONDS - B", 
                           "BONDS - Below B", 
                           "BONDS - Other")
      
        colnames(BLEND_returns_SMALL) <- lista_aktywow
        colnames(BLEND_returns_MEDIUM) <- lista_aktywow
        colnames(BLEND_returns_LARGE) <- lista_aktywow
        colnames(GROWTH_returns_SMALL) <- lista_aktywow
        colnames(GROWTH_returns_MEDIUM) <- lista_aktywow
        colnames(GROWTH_returns_LARGE) <- lista_aktywow
        colnames(VALUE_returns_SMALL) <- lista_aktywow
        colnames(VALUE_returns_MEDIUM) <- lista_aktywow
        colnames(VALUE_returns_LARGE) <- lista_aktywow
      
      
      # STEP 2 - Crating PORTFOLIO object
        Specyfikacja <- portfolio.spec(assets = lista_aktywow)
      
        
      # STEP 3 - Adding constraints constraints to portfolio
        # setting sum of weights in portfolio - in this case weight_sum = 1
          Specyfikacja <- add.constraint(portfolio=Specyfikacja, type="leverage", min_sum=0.9999, max_sum=1.0001)
        
        # Setting the weight sum for assets in portfolio - in this case 0,00 < asset_weight > 0,65
          Specyfikacja <- add.constraint(portfolio=Specyfikacja, type="box", min=0.0000, max=0.6500)
      
        
      # STEP 4 - Adding objectives in portfolio
        # Setting goal of maximalization of mean return
          CEL_Portfolio <- add.objective(portfolio=Specyfikacja, type='return', name='mean')  
        
        # Setting goal of minimalization of risk
          CEL_Portfolio <- add.objective(portfolio=CEL_Portfolio, type='risk', name='var')
      
        
      # STEP 5 - PORTFOLIO OPTIMIZATION
        # Portfolios - GROWTH Strategy
          Portfolio_GROWTH_SMALL <- optimize.portfolio(R = GROWTH_returns_SMALL, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          Portfolio_GROWTH_MEDIUM <- optimize.portfolio(R = GROWTH_returns_MEDIUM, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          Portfolio_GROWTH_LARGE <- optimize.portfolio(R = GROWTH_returns_LARGE, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          
        # Portfolios - VALUE Strategy
          Portfolio_VALUE_SMALL <- optimize.portfolio(R = VALUE_returns_SMALL, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          Portfolio_VALUE_MEDIUM <- optimize.portfolio(R = VALUE_returns_MEDIUM, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          Portfolio_VALUE_LARGE <- optimize.portfolio(R = VALUE_returns_LARGE, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          
        # Portfolios - BLEND Strategy
          Portfolio_BLEND_SMALL <- optimize.portfolio(R = BLEND_returns_SMALL, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          Portfolio_BLEND_MEDIUM <- optimize.portfolio(R = BLEND_returns_MEDIUM, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          Portfolio_BLEND_LARGE <- optimize.portfolio(R = BLEND_returns_LARGE, portfolio = CEL_Portfolio, optimize_method = "ROI", trace=TRUE)
          
        # Tabelaric summary of created portfolios with "ROI" optimization method
        
          Portfele_ROI_Tabela <- rbind(
                                       t(data.frame(Portfolio_GROWTH_SMALL[['weights']])),
                                       t(data.frame(Portfolio_GROWTH_MEDIUM[['weights']])),
                                       t(data.frame(Portfolio_GROWTH_LARGE[['weights']])),
                                      
                                       t(data.frame(Portfolio_VALUE_SMALL[['weights']])),
                                       t(data.frame(Portfolio_VALUE_MEDIUM[['weights']])),
                                       t(data.frame(Portfolio_VALUE_LARGE[['weights']])),
                                      
                                       t(data.frame(Portfolio_BLEND_SMALL[['weights']])),
                                       t(data.frame(Portfolio_BLEND_MEDIUM[['weights']])),
                                       t(data.frame(Portfolio_BLEND_LARGE[['weights']]))
                                      )
          
          Portfele_ROI_Tabela <- t(format(round(Portfele_ROI_Tabela, 4)))
          
          colnames(Portfele_ROI_Tabela) <- c(
                                             "GROWTH Strategy SMALL Capitalization",
                                             "GROWTH Strategy MEDIUM Capitalization",
                                             "GROWTH Strategy LARGE Capitalization",
                                             
                                             "VALUE Strategy SMALL Capitalization",
                                             "VALUE Strategy MEDIUM Capitalization",
                                             "VALUE Strategy LARGE Capitalization",
                                              
                                             "BLEND Strategy SMALL Capitalization",
                                             "BLEND Strategy MEDIUM Capitalization",
                                             "BLEND Strategy LARGE Capitalization"
                                            )
          
          write.table(Portfele_ROI_Tabela,"portfolios_ROI.csv",sep = "|")
    

# [06] *** =======[ BUSINESS CASES ]======= ***
        
  # Descriptions and outcomes
    print(Portfolio_GROWTH_SMALL)
    print(Portfolio_BLEND_MEDIUM)
    print(Portfolio_VALUE_LARGE)
    business_case_01 <- data.frame(Portfolio_GROWTH_SMALL[["weights"]])
    business_case_02 <- data.frame(Portfolio_BLEND_MEDIUM[["weights"]])
    business_case_03 <- data.frame(Portfolio_VALUE_LARGE[["weights"]])
    business_case_01 <- subset(business_case_01, business_case_01$Portfolio_GROWTH_SMALL...weights...>0.01)
    business_case_02 <- subset(business_case_02, business_case_02$Portfolio_BLEND_MEDIUM...weights...>0.01)
    business_case_03 <- subset(business_case_03, business_case_03$Portfolio_VALUE_LARGE...weights...>0.01)
  
    
  # Graphs for each Business Case
    
    # Business Case 01 - Individual Investor (GROWTH Strategy, SMALL Capitalization)
    
      png(file = "12 - Optimalization - GROWTH Strategy, SMALL Capitalization.jpg", width = 900, height = 900, units = "px")
      
      procenty_kolowe_2 <- sprintf("%.2f%s", round(business_case_01$Portfolio_GROWTH_SMALL...weights...*100 / sum(business_case_01$Portfolio_GROWTH_SMALL...weights...), 2), "%")
      pie(
          business_case_01$Portfolio_GROWTH_SMALL...weights..., 
          labels = procenty_kolowe_2, 
          main = "Optimalization - GROWTH Strategy, SMALL Capitalization",
          space = 0, 
          border = "white", 
          clockwise = TRUE, 
          col = c("darkslateblue", "yellow3")
      )
      
      symbols(0, 0, circles = 1, add=TRUE, bg="white", fg = "white", inches = 3)
      legend("topright",
             rownames(business_case_01),
             cex = 1.5,
             bty = "n",
             fill = c("darkslateblue", "yellow3"))
      
      dev.off()
    
      
    # Business Case 02 - Institution (BLEND Strategy, MEDIUM Capitalization)
      
      png(file = "13 - Optimalization - BLEND Strategy, MEDIUM Capitalization.jpg", width = 900, height = 900, units = "px")
      
      procenty_kolowe_3 <- sprintf("%.2f%s", round(business_case_02$Portfolio_BLEND_MEDIUM...weights...*100 / sum(business_case_02$Portfolio_BLEND_MEDIUM...weights...), 2), "%")
      pie(
          business_case_02$Portfolio_BLEND_MEDIUM...weights..., 
          labels = procenty_kolowe_3, 
          main = "Optimalization - BLEND Strategy, MEDIUM Capitalization",
          space = 0, 
          border = "white", 
          clockwise = TRUE, 
          col = c("seagreen","indianred")
        )
      
      symbols(0, 0, circles = 1, add=TRUE, bg="white", fg = "white", inches = 3)
      legend("topright",
             rownames(business_case_02),
             cex = 1.5,
             bty = "n",
             fill = c("seagreen","indianred"))
      
      dev.off()
    
      
    # Business Case 03 - Big corporation (VALUE Strategy, LARGE Capitalization)
      
      png(file = "14 - Optimalization - VALUE Strategy, LARGE Capitalization.jpg", width = 900, height = 900, units = "px")
      
      procenty_kolowe_4 <- sprintf("%.2f%s", round(business_case_03$Portfolio_VALUE_LARGE...weights...*100 / sum(business_case_03$Portfolio_VALUE_LARGE...weights...), 2), "%")
      pie(
          business_case_03$Portfolio_VALUE_LARGE...weights..., 
          labels = procenty_kolowe_4, 
          main = "Optimalization - VALUE Strategy, LARGE Capitalization",
          space = 0, 
          border = "white", 
          clockwise = TRUE, 
          col = c("peru","royalblue")
         )
      
      symbols(0, 0, circles = 1, add=TRUE, bg="white", fg = "white", inches = 3)
      legend("topright",
             rownames(business_case_03),
             cex = 1.5,
             bty = "n",
             fill = c("peru","royalblue"))
      dev.off()



























# WYKRESY GGPLOT2


# load library
library(ggplot2)

# Create test data.
data <- data.frame(
  category=c("A", "B", "C"),
  count=c(10, 60, 30)
)

# Compute percentages
data$fraction = data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n=-1))

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) # Try to remove that to see how to make a pie chart
























# INSPIRACJE DO MODELOWANIA


# 01 - Packages

install.packages("caret", dependencies = TRUE)
install.packager("mlbench")
install.packages("kernlab")
install.packages("randomForest")
install.packages("LiblineaR")
install.packages("ranger")
install.packages("naivebayes")
install.packages("kknn")
install.packages("HiDimDA")
install.packages("corrplot")

library(caret)
library(mlbench)
library(kernlab)
library(randomForest)
library(LiblineaR)
library(ranger)
library(naivebayes)
library(kknn)
library(HiDimDA)
library(corrplot)


# 02 - Data Set

data_set <- read.csv("Stroke_Prediction.csv", sep = ",", header = TRUE, na = c("NA", ""))

data_set <- na.omit(data_set)
data_set$Stroke <- as.factor(data_set$Stroke)


# 03 - Setting Data Set for Training and Testing

index_validation <- createDataPartition(data_set$Stroke, p=0.70, list=FALSE)

TRAIN_70 <- data_set[index_validation,]
TEST_30 <-data_set[-index_validation,]





# 05 - Models training and testing

control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
seed <- 7

# Linear Discriminant Analysis
set.seed(seed)
fit.lda <- train(Stroke~., data=TRAIN_70, method="lda", metric=metric, trControl=control)

# Decision Tree
set.seed(seed)
fit.cart <- train(Stroke~., data=TRAIN_70, method="rpart", metric=metric, trControl=control)

# kNN
set.seed(seed)
fit.knn <- train(Stroke~., data=TRAIN_70, method="knn", metric=metric, trControl=control)

# SVM - 1
set.seed(seed)
fit.svmLinearWeights <- train(Stroke~., data=TRAIN_70, method="svmLinearWeights", metric=metric, trControl=control)

# SVM - 2
set.seed(seed)
fit.svmLinear2 <- train(Stroke~., data=TRAIN_70, method="svmLinear2", metric=metric, trControl=control)

# Random Forest - 1
set.seed(seed)
fit.ranger <- train(Stroke~., data=TRAIN_70, method="ranger", metric=metric, trControl=control)

# Random Forest - 2
set.seed(seed)
fit.wsrf <- train(Stroke~., data=TRAIN_70, method="wsrf", metric=metric, trControl=control)

# Logistic Regression
set.seed(seed)
fit.regLogistic <- train(Stroke~., data=TRAIN_70, method="regLogistic", metric=metric, trControl=control)

# Naive Bayes
set.seed(seed)
fit.naive_bayes <- train(Stroke~., data=TRAIN_70, method = "naive_bayes", metric=metric, trControl=control)

# kkNN
set.seed(seed)
fit.kknn <- train(Stroke~., data=TRAIN_70, method="kknn", metric=metric, trControl=control)

# Factor-Based Linear Discriminant Analysis
set.seed(seed)
fit.RFlda <- train(Stroke~., data=TRAIN_70, method="RFlda", metric=metric, trControl=control)

# Neural network
set.seed(seed)
fit.nnet <- train(Stroke~., data=TRAIN_70, method="nnet", metric=metric, trControl=control)


# 06 - Accuracy summary

results <- resamples(list(lda = fit.lda, cart = fit.cart, knn = fit.knn,
                          svmLinearWeights = fit.svmLinearWeights, svmLinear2 = fit.svmLinear2, ranger = fit.ranger, wsrf = fit.wsrf,
                          regLogistic = fit.regLogistic, naive_bayes = fit.naive_bayes, kknn = fit.kknn, RFlda = fit.RFlda, nnet = fit.nnet))

summary(results)

# Accuracy comparison
dotplot(results)

# Best Model summary
fit.best = fit.svmLinear2
print(fit.best)

# Predictions and skill estimation of the best model on the testing dataset
predictions <- predict(fit.best, TEST_30)
confusionMatrix(predictions, TEST_30$Stroke)

