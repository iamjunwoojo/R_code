#산점도
ggplot(mtcars, aes(cyl, mpg)) +
  geom_point()
  
#점에 색깔 농도, size 입히기
ggplot(mtcars, aes(wt, mpg, size = disp, color = disp)) +
  geom_point()
  
#산점도 + smooth curve
ggplot(diamonds, aes(carat, price)) +
  geom_point() +
  geom_smooth()
  
#색깔별로 분류 color=clarity
ggplot(diamonds, aes(carat, price, color = clarity)) +
  geom_point() +
  geom_smooth()
  
#점 투명도 geom_point 에 있는 alpha 파라미터 
ggplot(diamonds, aes(carat, price, color = clarity)) +
  geom_point(alpha = 0.8) +
  geom_smooth()
  
  
#geom_point 의 color는 outline색깔만 지정 그래서 fill은 안에 색깔 but shape dfault에서 바꿔야됨
ggplot(mtcars, aes(wt, mpg, fill = fcyl)) +
  geom_point(shape = 21, size = 4, alpha = 0.6)

#geom_point의 aes 중 size로 분류
ggplot(mtcars,aes(mpg,wt,size=carb))+
  geom_point()
or
ggplot(mtcars,aes(mpg,wt))+
  geom_point(size=carb)
  
  
#alpha 에 변수넣으면 변수의 크기에 따른  투명도 지정가능  
plt_mpg_vs_wt <- ggplot(mtcars, aes(wt, mpg))


plt_mpg_vs_wt +
  geom_point(aes(alpha = fcyl))
  
#shape로 나누어 plot 만듬
ggplot(mtcars,aes(mpg,hp))+
  geom_point(aes(shape = fcyl))
  
#숫자로 라벨링
+geom_text(aes(label = fcyl))


#색깔지정
ggplot(mtcars, aes(wt, mpg)) +
  geom_point(color = my_blue, alpha = 0.6)
  
#텍스트레이어 붙이고 색깔지정
ggplot(mtcars, aes(wt, mpg, color = fcyl)) +
  geom_text(label = rownames(mtcars), color = "red")
  
#size 주목 저렇게도됨
ggplot(mtcars, aes(mpg, qsec, color = fcyl, shape = fam, size = hp/wt)) +
  geom_point()
  
  
#심플 bar plot
ggplot(mtcars, aes(fcyl, fill = fam)) +
  geom_bar() +
  labs(x = "Number of Cylinders", y = "Count")


#fill 색깔 지정
palette <- c(automatic = "#377EB8", manual = "#E41A1C")
ggplot(mtcars, aes(fcyl, fill = fam)) +
  geom_bar() +
  labs(x = "Number of Cylinders", y = "Count") +
  scale_fill_manual("Transmission", values = palette)
  
  
#fill legend 고치기
ggplot(mtcars,aes(gear,fill=fcyl))+
  geom_bar(position='dodge')+
  labs(x='gear',y='fcyl')+
  scale_fill_manual('fcyl 명단',values=c('yellow','green','blue'))

#zitter 이용해서 radom 하게 조금씩 변화시키고 y축 숫자로 지정

ggplot(mtcars, aes(mpg,0)) +
  geom_point(position = "jitter")

#zitter하고 y축 경계지정
  +ylim(c(-2, 2))


#zitter 하고 zitter 범위 지정 geom_jitter = position jitter +geom_point
plt_mpg_vs_fcyl_by_fam + geom_point()
plt_mpg_vs_fcyl_by_fam + geom_point(position = position_jitter(width = 0.1))
