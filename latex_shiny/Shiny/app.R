library(shiny)

ui <- fluidPage(
  
  titlePanel("Reparameterized in terms of the quantile Chen distribution"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      br(),
      
      sliderInput("lambda",
                  "λ",
                  value = 0.8,
                  min = 0.1,
                  max = 2,
                  step = 0.01),
      
      sliderInput("mu",
                  "μ",
                  value = 3,
                  min = 0.1,
                  max = 15,
                  step = 0.01),
      
      sliderInput("tau",
                  "τ",
                  value = 0.5,
                  min = 0.1,
                  max = 0.99,
                  step = 0.01)
      
    ),
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot")),
      )
    )
  )
)

server <- function(input, output) {
  d <- reactive({
    
    dist(input$lambda)
    dist(input$mu)
    dist(input$tau)
    
  })
  
  output$plot <- renderPlot({
    lambda <- input$lambda
    mu <- input$mu
    tau <- input$tau
    par(xaxs = "i", yaxs = "i")
    chen_reparametrizada <- function(y, lambda, mu, tau){
      exp1 <- (log(1-tau)/(1-exp(mu^lambda)))*lambda*(y^(lambda-1))
      exp2 <- exp((log(1-tau)/(1-exp(mu^lambda)))* (1-exp(y^lambda)) + y^lambda)
      return (exp1 * exp2)
    }
    
    curve(chen_reparametrizada(x, lambda, mu, tau) ,n=1000, from=0, to=20, ylim=c(0,0.5), ann=FALSE)
  })
}

shinyApp(ui, server)