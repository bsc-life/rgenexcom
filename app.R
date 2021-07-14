#########################################################################################
# SHINY APP THAT PERFORMS COMMUNITY DETECTION ALGORITHMS TO FIND META-PATIENT'S COMMUNITIES
# 
# V1
#
# Beatriz Urda García 2020
#########################################################################################

# shiny::runApp(system.file("shiny", package = "visNetwork"))
# setwd("Desktop/ANALYSIS/shiny_network_app/")

library(shiny) # runExample("01_hello")
library(igraph)
library(magrittr)
library(visNetwork)
library(data.table)
library(DT)
library(stringr)
library(shinydashboard)
library("shinydashboardPlus")
library(shinythemes)

# To share the library
library(rsconnect)
# rsconnect::deployApp('/home/beatriz/Desktop/ANALYSIS/shiny_network_app')

# net_filename <- 'final_pairwise_union_spearman_distance_sDEGs_network.txt'
#net_filename <- 'final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt'
# df <- read.csv(net_filename,header=T, sep="\t",stringsAsFactors = F); dim(df)
set.seed(5)

# df$lty = NA
# df[df$in_epidem == TRUE, ]$lty <- 1; df[df$in_epidem == FALSE, ]$lty <- 2

meta <- read.csv("new_disease_metadata_final_names.txt",header=T, sep="\t",stringsAsFactors = F); dim(meta)
discats <- unique(meta$disease_cat)

dis_colors <- meta[, c("final_disease_name","new_dis_cat_colors","disease_cat")]; colnames(dis_colors) = c("Dis", "Color","Dis_category")

discat_abbv <- data.frame(discats = sort(unique(meta$disease_cat)), 
                          abbv = c("congenital","circulatory","digestive","musculoskeletal",
                                   "nervous","respiratory","infectious","mentaldisorders","neoplasms"))

reactome_cats <- read.csv("aditional_data/reactome_pathway_categories.txt",header=F, sep="\t",stringsAsFactors = F)
reactome_cats <- reactome_cats$V1; dim(reactome_cats)
reactome_cats <- reactome_cats[!reactome_cats == "Digestion and absorption"]   # It is not present. 

diseases_dic <- readRDS("disease_dic.rds")
extended_node_names <- read.csv("aditional_data/extended_node_names.txt",header=T, sep="\t",stringsAsFactors = F); dim(extended_node_names)


# User Interface, generates the html
ui <- fluidPage(   # makes it responsible to different web browser windows dimensions
  # titlePanel("Community detection in disease similarity networks"),

  # tags$style(
  #   HTML('
  #       #legends {
  #       display: flex;
  #       align-items: center;
  #       justify-content: center;
  #       }
  #       #fluidrow1 {
  #       height:50px;
  #       }
  #      ')
  # ),
  
  navbarPage(title="From comorbidities to Gene Expression Fingerprints and back", #theme=shinytheme("flatly"),   # flatly, lumen
             tabPanel("Networks",   # NETWORKS SECTION
                      sidebarLayout(
                        sidebarPanel(
                          width = 2,
                          radioButtons("net_choice", "Select network:",
                                       choices=c("Disease similarity network (DSN)" = "final_pairwise_union_spearman_distance_sDEGs_network.txt",
                                                 "Stratified similarity network (SSN)" = "final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt")),
                          selectInput("pos_neg_interactions", "Select interactions:",
                                      choices=c("All" = "all",
                                                "Positive" = "pos",
                                                "Negative" = "neg")),
                          # radioButtons("pos_neg_interactions", "Select interactions:",
                          #              choices=c("Positive" = "pos",
                          #                        "Negative" = "neg")),
                          sliderInput("corr_slider","Edge's weight (|Spearman's correlation|):",min=0, max=1, value=c(0,1)),
                          radioButtons("comm_algorithm", "Color by:",
                                       choices=c("ICD9 category" = "ICD9",
                                                 "Greedy modularity optimization algorithm" = "greedy",
                                                 "Random walks" = "rand_walks")),  #walktrap
                          selectInput("filt_network_node", "Filter by disease(s)", sort(meta$final_disease_name), multiple=TRUE),
                          textOutput("net_choicePlot")
                        ),
                        mainPanel(
                          width = 10,
                          # textOutput("txt"),
                          # htmlOutput("txt"),
                          h2("Disease Similarity Networks based on gene expression profiles"),
                          p("Red edges denote posite interactions", style="color:#cc5052"),   # style="color:#cc5052; display:inline"
                          p("Blue edges denote negative interactions", style="color:#6b7ab5; margin-top: -0.4em"),
                          p("Dashed red lines denote positive interactions that correspond to known comorbidities (See Documentation)", style="margin-top: -0.4em"),
                          p("Nodes are coloured according to their ICD9 disease category", style="margin-top: -0.4em"),
                          br(),
                          fluidRow(id="fluidrow1",#div(style="display: inline-block;"),
                            column(width=9, offset = 0, aligh="center", visNetworkOutput("net_plot", height = "800px")),
                            column(width=3, offset = 0, align="center", id="legends",
                                   img(src="legend_discats_500.png", width='95%'))   # To center the image: , style="margin-top: +20em"
                          ),
                          # plotOutput("net_plot"),
                          br(), br(), br(),
                          # tableOutput("net_table"),
                          # fluidRow(
                          #   column(width=11, DTOutput("net_dt_table",width='80%'))),
                          
                          DTOutput("net_dt_table",width='95%'), # PRO Table
                          br(),
                          align="center"
                        )
                      )
             ),
             navbarMenu("Molecular mechanisms behind",    # MOL MECHANISMS SECTION
                        tabPanel("Diseases",
                                 sidebarLayout(
                                   sidebarPanel(
                                     width = 2,
                                     selectInput("reactome_cat_suppfig", "Select a Reactome pathway category: ", sort(reactome_cats))
                                   ),
                                   mainPanel(
                                     width = 10,
                                     align = 'center',
                                     h2("Reactome pathways significantly dysregulated in human diseases", align="center"),
                                     p("The heatmap shows the dysregulated Reactome pathways (rows) in the diseases (columns)"),   # style="color:#cc5052; display:inline"
                                     p("Overexpressed pathways", style="color:#6b7ab5; margin-top: -0.4em; display:inline"),
                                     p("and", style="margin-top: -0.4em; display:inline"),
                                     p("Underexpressed pathways", style="color:#cc5052; margin-top: -0.4em; display:inline"),
                                     br(),br(),
                                     fluidRow(align="center",
                                              imageOutput("fig_dis_fe", height = 'auto') #,width = 'auto'
                                     ),
                                     # br(),
                                     p("For each disease, Reactome pathways significantly up- and down-regulated were identified using the GSEA method (FDR <= 0.05)"),
                                     p("Ward2 algorithm was applied to cluster diseases based on the Euclidean distance of their binarized Normalized Effect Size (1s, and -1s for up- and down-regulated pathways)", 
                                       style="margin-top: -0.4em"),
                                     br(),br(),
                                   )
                                 )
                                 ),
                        tabPanel("Disease interactions",
                                 sidebarLayout(
                                   sidebarPanel(
                                     width = 2,
                                     h5("Select interactions between:"),
                                     selectInput("discats1", "Disease category: ",discats),
                                     selectInput("discats2", "Disease category: ",discats, selected = 'mental disorders')
                                   ),
                                   mainPanel(
                                     width = 10,
                                     h2("Over and Underexpressed pathways shared by", align="center"),
                                     h4("epidemiological interactions (EIs) and non-epidemiological interactions (NEIs)", align="center"),
                                     p("Each point represents a Reactome pathway category"),   # style="color:#cc5052; display:inline"
                                     p("The size of the points corresponds to the mean number of shared pathways in the EIs", style="margin-top: -0.4em"),
                                     p("The color corresponds to the ratio of the mean number of shared pathways in EIs versus NEIs", style="margin-top: -0.4em"),
                                     p("(e.g. red indicates that EIs share more altered pathways than NEIs)", style="margin-top: -0.4em"),
                                     br(),
                                     fluidRow(align="center",
                                              imageOutput("fig_dis_interactions") #,width = 'auto'
                                     ),
                                     br(),br(),
                                     align="center"
                                   )
                                 )
                        ),
                        tabPanel("Get genes and pathways",
                                 sidebarLayout(
                                   sidebarPanel(
                                     width = 2,
                                     # h5("Select a disease or diseases:"),
                                     selectInput("selected_node_genes", "Select disease(s):", sort(names(diseases_dic)), multiple=TRUE),
                                     selectInput("granularity", "Get: ",c("Genes", "Pathways")),
                                     radioButtons("granularity_signif", "Select features: ",
                                                  choices=c("All" = "all",
                                                            "Significantly altered" = "sign")),
                                     hr(),
                                     radioButtons("granularity_direction", NULL,
                                                  choices=c("All" = "all",
                                                            "Overexpressed" = "up",
                                                            "Underexpressed" = "down"))
                                   ),
                                   mainPanel(
                                     width = 10,
                                     align="center",
                                     h2("Genes and Pathways altered in diseases", align="center"),
                                     h4("or shared by disease groups and pairs", align="center"),
                                     br(),
                                     DTOutput("genes_pathways_table",width='95%'), # PRO Table
                                     DTOutput("common_genes_pathways",width='65%'), # PRO Table
                                     # textOutput("common_genes_pathways")
                                   )
                                 )
                        )
             ),
             tabPanel("Documentation",      # DOCUMENTATION SECTION
                      h2("From comorbidities to Gene Expression Fingerprints and back", align="left"),
                      h4("Beatriz Urda-García, Jon Sánchez-Valle, Rosalba Lepore and Alfonso Valencia", align="left"),
                      br(),
                      # column(width=3, offset = 0, aligh="left",
                      #        h4("Disease Similarity Network (DSN)", align="left"),
                      #        h5("First, we implemented an RNA-seq pipeline to perform Differential Expression Analysis for each disease.
                      #   Next, we computed the spearmans' correlation between the diseases' gene expression profiles. 
                      #   We kept significant interactions after multiple testing correction (FDR < 0.05)."),
                      #        h5("The obtained DSN contains positive and negative interactions, representing diseases with significantly 
                      # similar and dissimilar gene expression profiles, respectively.
                      #    Next, we evaluated the overlap of the positive interactions in the DSN with the epidemiological network from
                      #    Hidalgo et al. Positive interactions in the DSN described in this epidemiological network are represented with 
                      #    red dashed lines."),
                      #        br(),
                      #        h4("Stratified Similarity Network (SSN)", align="left"),
                      #        h5("We stratified each disease into subgroups of patients with similar expression profiles (meta-patients) by applying 
                      #    the k-medoids clustering algorithm to its normalized and batch effect corrected gene expression matrix.
                      #    Next, we performed Differential Expression Analyses for each meta-patient. "),
                      #        h5("To analyze the disease subtype-associated comorbidities, we built the Stratified Similarity Network (SSN) connecting 
                      #    meta-patients and diseases based on the similarities of their gene expression profiles 
                      #    (following the same methodology described for the DSN). The resulting Stratified Similarity Network (SSN) contains three 
                      #    types of interactions: (1) the previously described disease-disease interactions, (2) interactions connecting different 
                      #    meta-patients and (3) interactions connecting meta-patients to diseases."),
                      #        br()
                      #        ),
                      # column(width=3, offset = 0, aligh="right",
                      #        h4("Molecular mechanism behind diseases", align="left"),
                      #        h5("Reactome pathways significantly dysregulated in human diseases by pathway category.
                      # For each disease, Reactome pathways significantly over and underexpressed were identified using the GSEA method (FDR <= 0.05). 
                      # Ward2 algorithm was applied to cluster diseases based on the Euclidean distance of their binarized Normalized Effect Size 
                      # (1s, and -1s for over and underxpressed pathways). The heatmap shows the dysregulated Reactome pathways (rows) in the 
                      # diseases (columns), where over and underexpressed pathways are blue and red colored respectively."),
                      #        br(),
                      #        h4("Molecular mechanism behind disease interactions", align="left"),
                      #        h5("Over and underexpressed pathways behind epidemiological and not epidemiological interactions for each disease category pair.
                      #    Percentage of epidemiological versus non epidemiological interactions that share overexpressed or underexpressed 
                      #    pathways. Each point represents a Reactome pathway category. The size of the points corresponds to the mean number of shared
                      #    pathways in the epidemiological interactions. The color corresponds to the ratio of the mean number of shared pathways in 
                      #    epidemiological versus non epidemiological interactions."),
                      #        br(),
                      #        h4("Get genes and pathways", align="left"),
                      #        h5("In this section you can access the differentially expressed genes and pathways in a given phenotype (disease or meta-patient)
                      # and commonly dysregulated in phenotype pairs or groups."),
                      #        tags$li("If you select one phenotype, you will get the table of dysregulated genes and pathways
                      #    for that phenotype. You can filter the tables by selecting only the features that are significantly altered
                      #    or by selecting only the over or underexpressed features."),
                      #        tags$li("If you only select two or more phenotypes, you will get the genes or pathways
                      # that are significantly altered in all those phenotypes. Again, you can select only the over or the underexpressed features.")
                      # ),
                      
                      # h3("Networks", align="left"),
                      h4("Disease Similarity Network (DSN)", align="left"),
                      h5("First, we implemented an RNA-seq pipeline to perform Differential Expression Analysis for each disease.
                        Next, we computed the spearmans' correlation between the diseases' gene expression profiles. 
                        We kept significant interactions after multiple testing correction (FDR < 0.05)."),
                      h5("The obtained DSN contains positive and negative interactions, representing diseases with significantly 
                      similar and dissimilar gene expression profiles, respectively.
                         Next, we evaluated the overlap of the positive interactions in the DSN with the epidemiological network from
                         Hidalgo et al. Positive interactions in the DSN described in this epidemiological network are represented with 
                         red dashed lines."),
                      br(),
                      h4("Stratified Similarity Network (SSN)", align="left"),
                      h5("We stratified each disease into subgroups of patients with similar expression profiles (meta-patients) by applying 
                         the k-medoids clustering algorithm to its normalized and batch effect corrected gene expression matrix.
                         Next, we performed Differential Expression Analyses for each meta-patient. "),
                      h5("To analyze the disease subtype-associated comorbidities, we built the Stratified Similarity Network (SSN) connecting 
                         meta-patients and diseases based on the similarities of their gene expression profiles 
                         (following the same methodology described for the DSN). The resulting Stratified Similarity Network (SSN) contains three 
                         types of interactions: (1) the previously described disease-disease interactions, (2) interactions connecting different 
                         meta-patients and (3) interactions connecting meta-patients to diseases."),
                      br(),
                      # h3("Molecular mechanisms", align="left"),
                      h4("Molecular mechanism behind diseases", align="left"),
                      h5("It shows the Reactome pathways significantly dysregulated in human diseases by pathway category.
                      For each disease, Reactome pathways significantly over and underexpressed were identified using the GSEA method (FDR <= 0.05). 
                      Ward2 algorithm was applied to cluster diseases based on the Euclidean distance of their binarized Normalized Effect Size 
                      (1s, and -1s for over and underxpressed pathways). The heatmap shows the dysregulated Reactome pathways (rows) in the 
                      diseases (columns), where over and underexpressed pathways are blue and red colored respectively."),
                      br(),
                      h4("Molecular mechanism behind disease interactions", align="left"),
                      h5("Over and underexpressed pathways behind epidemiological and not epidemiological interactions for each disease category pair.
                         Percentage of epidemiological versus non epidemiological interactions that share overexpressed or underexpressed 
                         pathways. Each point represents a Reactome pathway category. The size of the points corresponds to the mean number of shared
                         pathways in the epidemiological interactions. The color corresponds to the ratio of the mean number of shared pathways in 
                         epidemiological versus non epidemiological interactions."),
                      br(),
                      h4("Get genes and pathways", align="left"),
                      h5("In this section you can access the differentially expressed genes and pathways in a given phenotype (disease or meta-patient)
                      and commonly dysregulated in phenotype pairs or groups."),
                      tags$li("If you select one phenotype, you will get the table of dysregulated genes and pathways
                         for that phenotype. You can filter the tables by selecting only the features that are significantly altered
                         or by selecting only the over or underexpressed features."),
                      tags$li("If you only select two or more phenotypes, you will get the genes or pathways
                      that are significantly altered in all those phenotypes. Again, you can select only the over or the underexpressed features.")
                      
                      ),
             tabPanel("Authors",
                      fluidRow(
                        align="center",
                        box(tags$a(imageOutput("beaimage"),href="https://www.bsc.es/es/urda-beatriz/publications",target="_blank"),
                            class="darkableImage",onmouseout="this.style.opacity=1;this.filters.alpha.opacity=100",
                            onmouseover="this.style.opacity=0.6;this.filters.alpha.opacity=60",
                            width = 2,align = "center",height = 2),
                        box(tags$a(imageOutput("jonimage"),href="https://www.bsc.es/es/sanchez-jon/publications",target="_blank"),
                            class="darkableImage",onmouseout="this.style.opacity=1;this.filters.alpha.opacity=100",
                            onmouseover="this.style.opacity=0.6;this.filters.alpha.opacity=60",
                            width = 2,align = "center",height = 2),
                        box(tags$a(imageOutput("albaimage"),href="https://orcid.org/0000-0002-9481-2557",target="_blank"),
                            class="darkableImage",onmouseout="this.style.opacity=1;this.filters.alpha.opacity=100",
                            onmouseover="this.style.opacity=0.6;this.filters.alpha.opacity=60",
                            width = 2,align = "center",height = 2),
                        box(tags$a(imageOutput("alfonsoimage"),href="https://www.icrea.cat/en/Web/ScientificStaff/alfonsovalencia-244256#researcher-nav",target="_blank"),
                            class="darkableImage",onmouseout="this.style.opacity=1;this.filters.alpha.opacity=100",
                            onmouseover="this.style.opacity=0.6;this.filters.alpha.opacity=60",
                            width = 2,align = "center",height = 2)
                      )
                    )),      # AUTHORS SECTION

)

server <- function(input, output) {
  
  # output$txt <- renderText({
  #   paste(input$current_node_id,"-",input$net_plot_selected)
  # })
  
  # output$txt <- renderUI({
  #   HTML(paste('<b>','Highlight node:','</b>'))
  #   # paste(input$current_node_id,"-",input$net_plot_selected)
  # })
  
  output$fig_dis_fe <- renderImage({
    cselection = tolower(gsub(" ", "_",input$reactome_cat_suppfig))
    cfilename = paste0("www/FE_dis/",cselection,".png")
    if(file.exists(cfilename)){
      list(src = cfilename, contentType = "image/png", width="70%", height="auto")
    }
  }, deleteFile = FALSE)
  
  output$fig_dis_interactions <- renderImage({
    cat1 <- discat_abbv[discat_abbv$discats == input$discats1, ]$abbv
    cat2 <- discat_abbv[discat_abbv$discats == input$discats2, ]$abbv
    cfilename1 = paste0("www/dis_interactions/",cat1,"_",cat2,".png")
    cfilename2 = paste0("www/dis_interactions/",cat2,"_",cat1,".png")
    if(file.exists(cfilename1)){
      list(src = cfilename1, contentType = "image/png", width="100%", height="auto")
    }else if(file.exists(cfilename2)){
      list(src = cfilename2, contentType = "image/png", width="100%", height="auto") #height=700
    }else{
      
    }
  }, deleteFile = FALSE)
  
  output$net_plot <- renderVisNetwork({
    df <- read.csv(input$net_choice,header=T, sep="\t",stringsAsFactors = F)
    if(input$pos_neg_interactions == 'pos'){
      # Selecting only positive interactions
      df <- df[df$Distance >= 0,]
    }else if(input$pos_neg_interactions == 'neg'){
      # Selecting only negative interactions
      df <- df[df$Distance < 0,]
      # df$Distance <- abs(df$Distance)
    }
    df <- df[((abs(df$Distance) >= input$corr_slider[1]) & (abs(df$Distance) <= input$corr_slider[2])), ]
    # p(input$net_plot_selected)
    if(!is.null(input$filt_network_node)){
      if(input$net_choice == "final_pairwise_union_spearman_distance_sDEGs_network.txt"){
        df <- df[(df$Dis1 %in% input$filt_network_node | df$Dis2 %in% input$filt_network_node),]
      }else{
        df <- df[(df$Corr_Dis1 %in% input$filt_network_node | df$Corr_Dis2 %in% input$filt_network_node),]
      }
    }
    
    graph <- graph_from_data_frame(df, directed=FALSE)
    # is_weighted(graph)
    E(graph)$weight <- abs(df$Distance)  # Uncomment maybe
    # E(graph)$lty <- df$lty

    # is_weighted(graph)
    #### Detect the communities and add the community as a vertex attribute
    if(input$comm_algorithm == 'ICD9'){
      # By ICD9 Code
      if(input$net_choice == "final_pairwise_union_spearman_distance_sDEGs_network.txt"){
        colors <- merge(data.frame(Dis=V(graph)$name), dis_colors, all.x = TRUE, all.y=FALSE)
        colors = colors[match(V(graph)$name, colors$Dis), ]
      }else{
        dis_order <- str_trim(gsub(" \\d+","",V(graph)$name))
        colors <- merge(data.frame(Dis=dis_order), dis_colors, all.x = TRUE, all.y=FALSE)
        colors = colors[match(dis_order, colors$Dis), ]
      }
      
      V(graph)$color <- as.character(colors$Color)
      V(graph)$Dis_category <- as.character(colors$Dis_category)
      nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, color = V(graph)$color, group=V(graph)$Dis_category)
    }else{
      # COMMUNITY DETECTION
      if(input$comm_algorithm == 'greedy'){
        # Using greedy optimization of modularity
        fc <- fastgreedy.community(graph)
        V(graph)$community <- fc$membership
      }else if(input$comm_algorithm == 'rand_walks'){
        # Using random walks
        fc <- cluster_walktrap(graph)
        V(graph)$community <- fc$membership #membership(fc)
      }
      # Visualize the communities
      nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, group = V(graph)$community)
      nodes <- nodes[order(nodes$id, decreasing = F),]
      
    }
    
    # Edge coloring and visualization
    edges <- get.data.frame(graph, what="edges")[1:2]
    # edges$color <- rep("lightgrey",length(edges$from))
    # Orginical red and blue
    edges$color <- df$Distance;
   
    edges$color[edges$color < 0] <- "#6b7ab595" #  "#405191" # BLUE
    edges$color[edges$color > 0] <- "#cc505295"   # "#C7535595" #"#d46e6e" # "#C75355" # "#C81E17" # RED
    
    edges$value <- abs(df$Distance)
    edges$dashes <- df$in_epidem; edges$dashes <- tolower(as.character(df$in_epidem))
    edges$dashes[edges$dashes == 'false'] <- "[6,0]"
    edges$dashes[edges$dashes == 'true'] <- "[6,15]"
    # visNetwork(nodes, edges)
    # edges$dashes = df$in_epidem   # dashes / shadow ; to invert: !df$in_epidem
    # edges$dashes = c(df$in_epidem, "[5,15]")
    # edges$dashes = c("[6,15]") # Todo dashed
    # edges$dashes = list("[6,15]",df$in_epidem) # No funciona
      # paste(df$in_epidem,"[5,15]",sep=",")
    
    # edges <- data.frame(from = c(1,2), to = c(1,3),dashes = c("[10,10,2,2]", "[2]"))
    
    visNetwork(nodes, edges) %>%
      visExport() %>%
      visOptions(highlightNearest = list(enabled=TRUE, degree=1,algorithm="hierarchical",labelOnly=FALSE), 
                 nodesIdSelection = list(enabled=TRUE, main="Highlight a node")) %>%   # list(enabled=TRUE, selected="BreastCancer")
      visIgraphLayout() %>%
      visInteraction(multiselect = T) %>%
      # visLegend() %>%
      # visLegend(position="right", main="Group") %>% # legend community detection
      # visEdges(shadow=list(color="#FFFAB5")) %>%
      visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_id', nodes.nodes);
            ;}")
  })
  # observeEvent(input$current_node_id, {
  #   visNetworkProxy("net_plot") %>%
  #     visGetSelectedNodes()
  #     # visGetNodes()
  # })
  observe({
    visNetworkProxy("net_plot") %>%
      visGetSelectedNodes()
    # visGetNodes()
  })
  
  # PRO TABLE
  output$net_dt_table <- renderDT({
    df <- read.csv(input$net_choice,header=T, sep="\t",stringsAsFactors = F)
    if(input$pos_neg_interactions == 'pos'){ # POS / NEG INTERACTIONS
      # Selecting only positive interactions
      df <- df[df$Distance >= 0,]
    }else if(input$pos_neg_interactions == 'neg'){
      # Selecting only negative interactions
      df <- df[df$Distance < 0,]
    }
    df <- df[((abs(df$Distance) >= input$corr_slider[1]) & (abs(df$Distance) <= input$corr_slider[2])), ]
    if(!is.null(input$filt_network_node)){
      if(input$net_choice == "final_pairwise_union_spearman_distance_sDEGs_network.txt"){
        df <- df[(df$Dis1 %in% input$filt_network_node | df$Dis2 %in% input$filt_network_node),]
      }else{
        df <- df[(df$Corr_Dis1 %in% input$filt_network_node | df$Corr_Dis2 %in% input$filt_network_node),]
      }
    }
    # if(!is.null(input$filt_network_node)){
    #   df <- df[(df$Dis1 %in% input$filt_network_node | df$Dis2 %in% input$filt_network_node),]
    # }
    df <- df[, c("Dis1","Dis2","Distance", "pvalue", "adj_pvalue", "in_epidem")]
    df$Distance <- round(df$Distance, 4)
    df$pvalue <- formatC(df$pvalue, format = "e", digits = 2); df$adj_pvalue <- formatC(df$adj_pvalue, format = "e", digits = 2)
    if(is.null(input$net_plot_selected) | input$net_plot_selected == ''){
      colnames(df) <- c("Disease 1", "Disease 2","Spearman's correlation", "p-value", "adj.p-value", "In epidemiology")
      df
    }else{
      filtered <- df[(df$Dis1 == input$net_plot_selected | df$Dis2 == input$net_plot_selected),]
      # filtered <- df[(df$Dis1 == input$networkid_selected | df$Dis2 == input$networkid_selected),]
      # df %>%
      # filter((Dis1 == input$dis_choice | Dis2 == input$dis_choice))
      colnames(filtered) <- c("Disease 1", "Disease 2","Spearman's correlation", "p-value", "adj.p-value", "In epidemiology")
      filtered
    }
    
  })
  
  # PRO TABLE: GENES OR PATHWAYS
  output$genes_pathways_table <- renderDT({
    if(length(input$selected_node_genes) == 1){
      
      old_dis_name = extended_node_names[extended_node_names$final_names == input$selected_node_genes, ]$old_names
      
      # old_dis_name = meta[meta$final_disease_name == input$selected_node_genes, ]$disease_name
      
      if(input$granularity == "Genes"){
        df <- read.csv(paste0("aditional_data/Genes/",old_dis_name,"_DEGs.txt"),header=T, sep="\t",stringsAsFactors = F)
        if(input$granularity_signif == "sign"){
          df <- df[df$adj.p.val < 0.05, ]
        }else{df}
        if(input$granularity_direction == 'up'){
          df <- df[df$logFC > 0, ]
        }else if(input$granularity_direction == 'down'){
          df <- df[df$logFC < 0, ]
        }else{df}
        
        df$adj.p.value = formatC(df$adj.p.value, format = "e", digits = 2)
        df$logFC = round(df$logFC,3)
        df
        
      }else{
        df <- read.csv(paste0("aditional_data/FE/",old_dis_name,"_pathways.txt"),header=T, sep="\t",stringsAsFactors = F)
        if(input$granularity_signif == "sign"){
          df <- df[df$FDR.q.val < 0.05, ]
        }else{df}
        if(input$granularity_direction == 'up'){
          df <- df[df$NES > 0, ]
        }else if(input$granularity_direction == 'down'){
          df <- df[df$NES < 0, ]
        }else{df}
        df$NES = round(df$NES, 3)
        df$FDR.q.val = formatC(df$FDR.q.val, format = "e", digits = 2)
        df
      }
      
    }else if(length(input$selected_node_genes) > 1){
      
    }
  })
  
  output$common_genes_pathways <- renderDT({
    nsel_nodes = length(input$selected_node_genes)
    if(nsel_nodes > 1){
      intersection <- c()
      if(input$granularity == "Genes"){  # GENES
        
        if(input$granularity_direction == 'up'){ # up
          for(dis in input$selected_node_genes){
            if(dis == input$selected_node_genes[1]){
              intersection <- diseases_dic[[dis]]@genes_up
            }else{
              intersection = intersect(intersection, diseases_dic[[dis]]@genes_up)
            }
          }
          
        }else if(input$granularity_direction == 'down'){     # down
          for(dis in input$selected_node_genes){
            if(dis == input$selected_node_genes[1]){
              intersection <- diseases_dic[[dis]]@genes_down
            }else{
              intersection = intersect(intersection, diseases_dic[[dis]]@genes_down)
            }
          }
          
        }else{
          for(dis in input$selected_node_genes){
            if(dis == input$selected_node_genes[1]){
              intersection1 <- diseases_dic[[dis]]@genes_up
              intersection2 <- diseases_dic[[dis]]@genes_down
            }else{
              intersection1 = intersect(intersection1, diseases_dic[[dis]]@genes_up)
              intersection2 = intersect(intersection2, diseases_dic[[dis]]@genes_down)
            }
          }
          intersection = union(intersection1, intersection2)
          
        }
        
      }else{  # PATHWAYS
        
        if(input$granularity_direction == 'up'){   # up
          
          for(dis in input$selected_node_genes){
            if(dis == input$selected_node_genes[1]){
              intersection <- diseases_dic[[dis]]@pathways_up
            }else{
              intersection = intersect(intersection, diseases_dic[[dis]]@pathways_up)
            }
          }
          
        }else if(input$granularity_direction == 'down'){   # down
          for(dis in input$selected_node_genes){   
            if(dis == input$selected_node_genes[1]){
              intersection <- diseases_dic[[dis]]@pathways_down
            }else{
              intersection = intersect(intersection, diseases_dic[[dis]]@pathways_down)
            }
          }
          
        }else{
          for(dis in input$selected_node_genes){   
            if(dis == input$selected_node_genes[1]){
              intersection1 <- diseases_dic[[dis]]@pathways_up
              intersection2 <- diseases_dic[[dis]]@pathways_down
            }else{
              intersection1 = intersect(intersection1, diseases_dic[[dis]]@pathways_up)
              intersection2 = intersect(intersection2, diseases_dic[[dis]]@pathways_down)
            }
          }
          intersection = union(intersection1, intersection2)
        }
        
        
      }
      nintersect = length(intersection)
      data.frame(intersection)
      # paste0("Common significantly dysregulated genes: ", as.character(length(intersection)))
    }
  })
  
  #### AUTHORS IMAGES
  output$beaimage <- renderImage({
    list(src = "www/beaimage.png", contentType = "image/png", width="90%", height="auto")
  }, deleteFile = FALSE)
  output$jonimage <- renderImage({
    list(src = "www/jonimage.png", contentType = "image/png", width="90%", height="auto")
  }, deleteFile = FALSE)
  output$albaimage <- renderImage({
    list(src = "www/albaimage.png", contentType = "image/png", width="90%", height="auto")
  }, deleteFile = FALSE)
  output$alfonsoimage <- renderImage({
    list(src = "www/alfonsoimage.png", contentType = "image/png", width="90%", height="auto")
  }, deleteFile = FALSE)
  
}


shinyApp(ui = ui, server = server)


