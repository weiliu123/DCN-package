getDiffSubnetwork <-
function(graph, seed, d = 2, r = 0.05, alpha = 0.7, th = 0.5^6){
  adjm <- as_adjacency_matrix(graph, type = "both", attr = "DC", sparse = FALSE)  # 邻接矩阵
  v.DE <- vertex_attr(graph)$DE
  names(v.DE) <- vertex_attr(graph)$name
  subg <- make_empty_graph(0,directed = FALSE) %>% add_vertices(1, name = seed, DE = v.DE[seed])
  cur.v <- get.vertex.attribute(subg,"name")  #当前子图中的节点
  if(max(adjm[cur.v,]) > th){
    max.v <- colnames(adjm)[which.max(adjm[cur.v,])]
    subg <- add.vertices(subg, nv=1, name=max.v, DE = v.DE[max.v])
    subg <- add.edges(subg, edges = c(cur.v, max.v), attr = list(DC = adjm[cur.v, max.v]))
    adjm[cur.v, max.v] <- 0
    adjm[max.v, cur.v] <- 0
    
    DE <- getDE(subg)
    DC <- getDC(subg)
    D <- alpha*DE + (1-alpha)*DC
    
    cur.v <- max.v
    repeat{
      # 加当前网络中顶点之间的边
      # cur.v <- max.v
      v.subg <- get.vertex.attribute(subg,"name")   # 获取当前网络顶点名
      cur.adjm <- adjm[v.subg, v.subg]           # 去除subg.temp图中边之后的当前顶点的差异共表达网络邻接矩阵
      cur.v.adjv <- colnames(cur.adjm)[which(cur.adjm[cur.v,] > th)]
      if(length(cur.v.adjv) > 0){
        for(i in 1:length(cur.v.adjv)){
          subg <- add.edges(subg, edges = c(cur.v, cur.v.adjv[i]), attr = list(DC = cur.adjm[cur.v, cur.v.adjv[i]]))
        }
        DE <- getDE(subg)
        DC <- getDC(subg)
        D <- alpha*DE + (1-alpha)*DC
      }			
      
      # 是否需要把subg中顶点的邻接矩阵赋值为0
      adjm[v.subg, v.subg] <- 0
      
      # 增加与当前网路邻接的顶点
      
      # 计算各节点到种子节点距离
      distToSeed <- distances(subg, to=seed, weights=NA, algorithm="unweighted")
      
      # 选择到种子节点距离 < d 的节点作为候选节点，以增加新的邻接顶点
      v.candidate <- rownames(distToSeed)[which(distToSeed < d)]                  
      subg.adjm <- adjm[v.candidate, ]
      flag <- 0      # 记录是否加入新顶点 
      repeat{
        if(max(subg.adjm) > th){                     # 差异共表达网络阈值
          # 找与当前网络中顶点相连的具有最大权值的边，以及对应的顶点
          max.ind <- which(subg.adjm==max(subg.adjm), arr.ind = TRUE)     
          # browser()
          v1 <- v.candidate[max.ind[1,1]]                                      # 当前网络中的顶点
          v2 <- colnames(subg.adjm)[max.ind[1,2]]                         # 新顶点
          
          # 增加顶点
          subg.temp <- add.vertices(subg, nv=1, name = v2, DE = v.DE[v2])
          
          # 增加边
          subg.temp <- add.edges(subg.temp, edges = c(v1, v2), attr = list(DC = adjm[v1, v2]))  
          
          subg.adjm[v1, v2] <- 0
          
          cur.DE <- getDE(subg.temp)
          cur.DC <- getDC(subg.temp)
          cur.D <- alpha*cur.DE + (1-alpha)*cur.DC
          # browser()
          if(cur.D > (1+r)*D){
            D <- cur.D              # 更新DC
            DE <- cur.DE
            DC <- cur.DC					
            subg <- subg.temp         # 更新subg
            adjm[v1, v2] <- 0      # 更新邻接矩阵
            adjm[v2, v1] <- 0
            cur.v <- v2               # 更新cur.v, 以便下一步循环加入与该顶点相连的边
            # browser()
            flag <- 1    # 加入新顶点标志赋值为1
            break
          }
        }else{
          break
        }
      }
      if(flag == 0){
        break
      }
    }
  }
  
  
  DE <- getDE(subg)
  DC <- getDC(subg)
  D <- alpha*DE + (1-alpha)*DC
  
  return(list(subg = subg, D = D, DE = DE, DC = DC))
}
