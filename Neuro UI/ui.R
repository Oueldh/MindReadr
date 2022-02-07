
library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Visualizing EEG Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            
            #The file type
            radioButtons("data_type","Type of Data",c(csv=',',txt="")),
            
            #The Header in your data if it is present or not 
            radioButtons("h","Header",c(yes=T,no=F)),
            
            #Input files 
            fileInput("EEG_Closed","Eyes Closed"),
            fileInput("EEG_Open","Eyes Open"),
            
            #Threshold for the adjacency matrix
            sliderInput("thresh",
                        "Threshold:",
                        min = 0,
                        max = 1,
                        value =.5),
            #Sliding window
            numericInput("interval","Time Block",value = 5),
            
            #name of the video:
            textInput("name","Video name",value=""),
            #The band that the user wants
            selectInput("band","Band",choices = list("delta"=1,"theta"=2,
                                                    "alpha"=3,"beta"=4,
                                                    "gamma"=5)),
            
            
            #Submission button
            actionButton("go","Submit",class="btn btn-primary")
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           uiOutput("distPlot")
        )
    )
))
