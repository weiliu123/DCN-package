getDiffSubnetwork <-
function(graph, seed, d = 2, r = 0.05, alpha = 0.7, th = 0.5^6){
  adjm <- as_adjacency_matrix(graph, type = "both", attr = "DC", sparse = FALSE)  # �ڽӾ���
  v.DE <- vertex_attr(graph)$DE
  names(v.DE) <- vertex_attr(graph)$name
  subg <- make_empty_graph(0,directed = FALSE) %>% add_vertices(1, name = seed, DE = v.DE[seed])
  cur.v <- get.vertex.attribute(subg,"name")  #��ǰ��ͼ�еĽڵ�
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
      # �ӵ�ǰ�����ж���֮��ı�
      # cur.v <- max.v
      v.subg <- get.vertex.attribute(subg,"name")   # ��ȡ��ǰ���綥����
      cur.adjm <- adjm[v.subg, v.subg]           # ȥ��subg.tempͼ�б�֮��ĵ�ǰ����Ĳ��칲���������ڽӾ���
      cur.v.adjv <- colnames(cur.adjm)[which(cur.adjm[cur.v,] > th)]
      if(length(cur.v.adjv) > 0){
        for(i in 1:length(cur.v.adjv)){
          subg <- add.edges(subg, edges = c(cur.v, cur.v.adjv[i]), attr = list(DC = cur.adjm[cur.v, cur.v.adjv[i]]))
        }
        DE <- getDE(subg)
        DC <- getDC(subg)
        D <- alpha*DE + (1-alpha)*DC
      }			
      
      # �Ƿ���Ҫ��subg�ж�����ڽӾ���ֵΪ0
      adjm[v.subg, v.subg] <- 0
      
      # �����뵱ǰ��·�ڽӵĶ���
      
      # ������ڵ㵽���ӽڵ����
      distToSeed <- distances(subg, to=seed, weights=NA, algorithm="unweighted")
      
      # ѡ�����ӽڵ���� < d �Ľڵ���Ϊ��ѡ�ڵ㣬�������µ��ڽӶ���
      v.candidate <- rownames(distToSeed)[which(distToSeed < d)]                  
      subg.adjm <- adjm[v.candidate, ]
      flag <- 0      # ��¼�Ƿ�����¶��� 
      repeat{
        if(max(subg.adjm) > th){                     # ���칲����������ֵ
          # ���뵱ǰ�����ж��������ľ������Ȩֵ�ıߣ��Լ���Ӧ�Ķ���
          max.ind <- which(subg.adjm==max(subg.adjm), arr.ind = TRUE)     
          # browser()
          v1 <- v.candidate[max.ind[1,1]]                                      # ��ǰ�����еĶ���
          v2 <- colnames(subg.adjm)[max.ind[1,2]]                         # �¶���
          
          # ���Ӷ���
          subg.temp <- add.vertices(subg, nv=1, name = v2, DE = v.DE[v2])
          
          # ���ӱ�
          subg.temp <- add.edges(subg.temp, edges = c(v1, v2), attr = list(DC = adjm[v1, v2]))  
          
          subg.adjm[v1, v2] <- 0
          
          cur.DE <- getDE(subg.temp)
          cur.DC <- getDC(subg.temp)
          cur.D <- alpha*cur.DE + (1-alpha)*cur.DC
          # browser()
          if(cur.D > (1+r)*D){
            D <- cur.D              # ����DC
            DE <- cur.DE
            DC <- cur.DC					
            subg <- subg.temp         # ����subg
            adjm[v1, v2] <- 0      # �����ڽӾ���
            adjm[v2, v1] <- 0
            cur.v <- v2               # ����cur.v, �Ա���һ��ѭ��������ö��������ı�
            # browser()
            flag <- 1    # �����¶����־��ֵΪ1
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