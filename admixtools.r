#ADMIXTOOLS 比DSUITE 准确
#vcf文件生产map和ped文件
plink --vcf cyr34.LD.vcf --make-bed --allow-extra-chr --recode --autosome-num 18 --double-id --out test
#染色体名称需要替换为数字
sed -i "s/CM0458//g" test.map
sed -i "s/JANSJV0100000//g" test.map

#vim het_maf.parfile
genotypename: test.ped
snpname:      test.map
indivname:    test.ped
outputformat:    EIGENSTRAT
genooutfilename:   test.geno
snpoutfilename:    test.snp
indoutfilename:    test.ind
familnames: NO
#convertf转换格式
convertf -p het_maf.parfile

#R设置分组
Outgroup = 'Outgroup'
Europe = c('Europe_1', 'Europe_2', 'Europe_3', 'Europe_4', 'Europe_5')
Asia = c('Asia_1', 'Asia_2')
Africa = c('Central_Africa', 'South_Africa')
North_America = 'North_America'
South_America = 'South_America'
Oceania = 'Oceania'
pop1 <- c('Asia_1', 'Asia_2','Central_Africa', 'South_Africa','North_America','South_America','Oceania')
pop2 <- c('Asia_1', 'Asia_2','Central_Africa', 'South_Africa','North_America','South_America','Oceania')


#R导入数据，只放.geno .ind .snp三个文件，其中.ind中的？？？替换为分组信息
#auto_only = FALSE 当snp数目太少，加此参数，否则无法运行
f2_blocks<-f2_from_geno("data/test2",auto_only = FALSE)
f3 <- f3(f2_blocks, Europe, pop1, pop2)
f4 <- f4(f2_blocks, pop1, pop2, Europe, Outgroup, f4mode = FALSE)

#输出文件
write.table(f3,file = "f3.txt",quote=F,col.name=T,sep ="\t",row.names=F)
write.table(f4,file = "f4.txt",quote=F,col.name=T,sep ="\t",row.names=F)

#f3的值应大部分大于0
#f3(A;B,C)=f4(A,B;A,C)=1/2(f2(A,B)+f2(A,C)−f2(B,C))
#根据公式可以看出f3(A;B,C) 可以检测A是否为B,C混合，A为B,C混合时，f3<0。
#也可以把A固定，比如A设置为外群来比较B,C的关系。由于A为外群，f2(A,B)，f2(A,C)基本是一致的，B和C种群越密切，f2(B,C))越小，f3(A;B,C)越大。

#f4结果解读（>0或者<0，结果相似）。
#D>0，P1和P3之间有基因流；D<0,P2和P3之间有基因流。

#计算fst
fst <- fst(f2_blocks, pop1 = c('Asia_1', 'Asia_2','Central_Africa', 'South_Africa','North_America','South_America','Oceania','Europe_1', 'Europe_2', 'Europe_3', 'Europe_4', 'Europe_5') , 
           pop2 = c('Asia_1', 'Asia_2','Central_Africa', 'South_Africa','North_America','South_America','Oceania','Europe_1', 'Europe_2', 'Europe_3', 'Europe_4', 'Europe_5'),
           adjust_pseudohaploid = FALSE)
write.table(fst,file = "fst.txt",quote=F,col.name=T,sep ="\t",row.names=F)
