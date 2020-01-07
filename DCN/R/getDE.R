getDE <-
function(g){
  v.DE <- get.vertex.attribute(g,"DE")
  DE <- sum(abs(v.DE)) / sqrt(length(v.DE))  #计算差异表达得分
}
