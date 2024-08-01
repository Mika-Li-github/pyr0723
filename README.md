使用 biomaRt 存取 Ensembl 註釋
http://127.0.0.1:26012/library/biomaRt/doc/accessing_ensembl.html

20240730:
今天果然還是先不隱藏變數。
流程:
把網頁資料(X2)整合
pandasGWAS再抓一次，再刪除不在基因上的
其他不夠的先用biomaRt找
真的找不到就去TW biobank

在ensemble發現了Phenotype/Disease/Trait這一欄，之後找大概有哪些，再擴大搜尋
把GWAS查到的Traits，想要的都存下來了
做好兩個bigPanda，現在要去R合併

20240801:
先手動刪除載體資料
Asso沒有chr資料
