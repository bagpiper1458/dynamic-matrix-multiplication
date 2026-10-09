# 從原始公式推導 PARAMS 的識別、block、預計算與生命週期

本文件從作者提供的論文 Algorithms 1–4、一般公式，以及後來找到的 lowbit 公式產生器推導一套可以實作的 PARAMS 語義。這是對 pseudo code 未展開細節的具體化，不宣稱是作者唯一可能的原意，也不宣稱全域最優。

此階段沒有更改既有 Python/C++ 運算核心或 benchmark。新增的 `tools/derive_params.py` 是小尺寸符號／數值核對工具，不是要取代效能版本。

## 共同記號：先定義「通式的項」

令 N=2ⁿ，i 是輸出列、j 是通式的項遮罩，範圍都是 0…N−1；補集只取 n 個 bits。GF(2⁸) 加法為 XOR。

\[
X(i,j)=\bigoplus_{s\subseteq j} a_{(i\mathbin{\&}j)-s},
\qquad
Y(i,j)=\bigoplus_{t\subseteq((N-1)\oplus j)}b_{(j\oplus(i\mathbin{\&}j))\,|\,t}.
\]

X 是原稿 SM_A 的集合和，Y 是 SM_B **完成欄反轉之後**的集合和。這與舊碼 `B[i][0]=N-1-i` 並倒序讀取 B 一致。circulant 將係數 a 的索引取 mod N。

定義「通式項表」M：

\[
M[i,j]=X(i,j)Y(i,j),\qquad z_i=\bigoplus_{j=0}^{N-1}M[i,j].
\]

M 不是原始矩陣 A；它是一般公式中的乘積排成的表。實作不必真的配置這張 N×N 表。

## 1. pos：識別的是子式，不是目前數值或出現座標

### 向量和的 pos

設 F 為自由位元遮罩，r 為固定為 1 的位元，要求 F&r=0，且兩者都不超過 n bits：

\[
V(F,r)=\bigoplus_{s\subseteq F}b_{r|s},\qquad
\operatorname{pos}_B=(\text{B},F,r).
\]

對通式的 Y(i,j)，F=(N−1) XOR j，r=j XOR (i&j)。例如 N=4：V(01₂,00₂)=b₀⊕b₁，V(11₂,00₂)=b₀⊕b₁⊕b₂⊕b₃。

每個 bit 有「自由、固定 0、固定 1」三種狀態，所以不同向量和共 3ⁿ 個，其中 N 個是原始 b 元素。

計算時選 F 中的一個位元 q（本排程固定選最高位元）：

\[
V(F,r)=V(F\setminus\{q\},r)\oplus V(F\setminus\{q\},r|q).
\]

每個非原始元素只需計算一次，因此共 3ⁿ−N 次 XOR。所有 dependency 先完成，再使用其 value。

### 係數與乘積的 pos

一般 Toeplitz 可將係數和表示為各 a 索引出現次數的奇偶 bitset；相同索引出現兩次就取消。circulant 先把索引取 mod N，再做奇偶正規化。這是精確的符號識別方式，不是比較計算後 byte 值是否相同。

```text
pos_A = (A, normalized_parity_of_coefficient_indices)
pos_P = (P, pos_A, pos_B)
```

這是容易核對的參考表示。正式執行可以將 pos 編譯成緊密整數 ID，避免每個向量都重建 bitset 或查 hash。現有係數 DAG 的 canonical operand-pair 去重是一種較便宜、但不完備的符號共用方式；它不會識別所有可能的代數等價式。

不同子式恰好算出同一個 byte，不代表可以在下一個輸入沿用。不同尺寸／矩陣種類的 layout 也不可混用；以下 pos 均以固定的 layout 為作用域。

## 2. block：對齊的方塊，value 是各列部分和

選定以原點對齊、邊長 h=2ᵏ 的區塊。R、C 是起始列／欄，均為 h 的倍數。定義：

\[
W_k(R,C)[u]=\bigoplus_{v=0}^{h-1}M[R+u,C+v],\qquad 0\le u<h.
\]

所以 1×1 block 的 value 是一個乘積，2×2 block 的 value 是兩個 byte，完整 N×N block 的 value 就是 N 個輸出。若使用連續陣列，value 可實作為一段 span，不必每個 block 都配置 vector。

這是對論文「處理 2ᵏ×2ᵏ block」的具體定義。若把整個 block 壓成單一 XOR，會丟失各輸出列，因此不能採用那個解讀。

### 不用展開 block，也能直接得到 pos

\[
m=C/2^k,\quad q=(R/2^k)\mathbin{\&}m,
\qquad\boxed{\operatorname{pos}_W=(\text{W},k,m,q)}.
\]

理由：將 i=R+u、j=C+v 代入 X/Y，低 k bits 由區塊內的 u、v 決定；高 bits 只依賴 m 及 q。m bit=0 時，對應的 R bit 不影響 X 或 Y，因此不同座標會落到同一個 pos。

每個高 bit 有 (m,q)=(0,0)、(1,0)、(1,1) 三種狀態；第 k 層因此只有 **3ⁿ⁻ᵏ 個結構 pos**。沒有必要走遍全部 N² 個出現位置才去重。

### 父區塊只依賴三個不同子區塊

令子區塊邊長 h/2：

```text
F = W[k-1](R,       C)
G = W[k-1](R,       C+h/2)
H = W[k-1](R+h/2,   C+h/2)
```

左下 W[k−1](R+h/2,C) 的 pos 與 F 完全相同，因為其 column mask 的新增 bit=0。因此：

\[
W_k(R,C)=\operatorname{concat}(F\oplus G,\ F\oplus H).
\]

對父 pos=(k,m,q)，三個子 pos 直接是：

```text
F: (k-1, 2m,   2q)
G: (k-1, 2m+1, 2q)
H: (k-1, 2m+1, 2q+1)
```

這是**通式項表的 block 關係**。它與論文對原始矩陣所做的 F/G/H 分解具有相同的合併結構；兩者的 F/G/H 記號在此表達對應子結果，而非宣稱項表等於原矩陣。

若以三進位整數 t 編碼每個 (m,q) 的高位元狀態，三個孩子恰好是 3t、3t+1、3t+2。這就推導出目前迭代版三個相鄰區塊的合併順序。

### 4×4 的直接例子

將不同乘積依現有三進位順序命名為 p₀…p₈，通式項表為：

\[
M=\begin{bmatrix}
p_0&p_1&p_3&p_4\\
p_0&p_2&p_3&p_5\\
p_0&p_1&p_6&p_7\\
p_0&p_2&p_6&p_8
\end{bmatrix}.
\]

四個對齊 2×2 區塊中，左上和左下是同一個 F。只需算：

```text
F = [p0 XOR p1, p0 XOR p2]     pos=(1,0,0)
G = [p3 XOR p4, p3 XOR p5]     pos=(1,1,0)
H = [p6 XOR p7, p6 XOR p8]     pos=(1,1,1)
result = concat(F XOR G, F XOR H)
```

不同乘積是 9 個；此層用 6 XOR，最後合併用 4 XOR，加上向量和的 5 XOR，得到論文的 **9M、15A**。係數準備另外計算。

### 一般的運算次數

第 k 層有 3ⁿ⁻ᵏ 個 block，每個產生長度 2ᵏ 的向量、需要 2ᵏ 次 XOR：

\[
\sum_{k=1}^{n}2^k3^{n-k}=2(3^n-N).
\]

加上向量和的 3ⁿ−N 次 XOR，正好是 3(3ⁿ−N)。乘積層是 3ⁿ 次 GF 乘法。unit-sum circulant 的 j=N−1 對應 N 個係數全為 1 的乘積，可直接複製向量元素，乘法降為 3ⁿ−N。

## 3. 預計算：分成 layout、矩陣、每個向量三種生命週期

| 資料 | 依賴什麼 | 何時可重用／何時失效 |
|---|---|---|
| C/L、索引、pos→ID、block 三子關係、讀寫順序 | N 與矩陣結構 | 同尺寸／結構可共用 |
| 每個值的最後使用點、槽位配置 | 已選定的排程 | layout 不變即可共用 |
| X 的係數值、共用係數子式 | 矩陣 A | A 改變就更新 |
| unit-sum 驗證結果及可省略乘法的位置 | 結構與 A 的有效條件 | A 改變需重新驗證 |
| Y 向量和 | b | b 改變就更新 |
| 乘積與 block 各列部分和 | A 與 b | 任一改變就更新 |

因此 `PARAMS` 的**位置、相依關係與排程**可以跨輸入保留；其中的 **value** 不一定可以。共用同一矩陣的多個向量時，保留係數值即可，不能把上一個向量的乘積當作有效快取。

若 value 用 byte 儲存，0…255 都是合法值，不能靠 byte 的 -1/255 判斷未初始化。可使用獨立 valid flag、epoch，或讓拓樸順序保證讀之前一定已寫入。原稿 value=-1 若採較寬的整數型別則可作 sentinel，但這是型別設計選擇。

## 4. 何時釋放／覆寫：按相依關係，而不是按插入時間

一般規則是：某值所有尚未完成的讀取都結束後，才能回收它的槽位。layout 編譯時可以預計算 last-use，或計算各節點剩餘讀取次數。F 同時供上下半部使用，必須讀完兩處才可回收。

例如 4×4 中 p₀ 同時供 F 的兩個分量使用，不能算完 p₀⊕p₁ 就把 p₀ 覆寫；需要先保存或讀完另一個 p₀⊕p₂。資料的符號 ID 永遠不變，但其 runtime 槽位在死亡後可以交給另一個 ID。

### 能直接實作的逐階段排程

1. **向量和展開**：依固定順序產生全部 3ⁿ 個 V。使用反向掃描，把 [u,v] 展成 [u⊕v,v,u]，避免覆寫尚未讀取的輸入。原始輸入不需要修改，先複製到 work。
2. **乘積**：所有 V 的相互依賴已計算完，每個 V 只剩對應乘積需要它，因此可直接 `work[t] = X[t] * work[t]`。unit 葉保留該值即可。
3. **block 合併**：在第 k 層，每個父節點讀三段長度 h=2ᵏ⁻¹ 的 F/G/H，輸出兩段。由前往後讀，每個 offset 的 f/g/h 都先載入，再寫 f⊕g 與 f⊕h。

合併一層的安全 in-place 索引為：

```text
source = group * 3h
target = group * 2h
for offset = 0 .. h-1:
    f = work[source + offset]
    g = work[source + h + offset]
    hvalue = work[source + 2h + offset]
    work[target + offset]     = f XOR g
    work[target + h + offset] = f XOR hvalue
```

group、offset 都遞增，三個來源先載入，才進行兩次寫入。這也是現有 iterative 的 contraction 規則；這次推導補上它與論文 blocks 的精確索引關係。

### 空間可以證明到哪裡

第 k 層所有唯一 block 的有效值合計：

\[
S_k=2^k3^{n-k},\qquad S_{k+1}=\tfrac23S_k.
\]

如果把所有層都永久存在 PARAMS，光 block 的 byte values 就有：

\[
\sum_{k=0}^nS_k=3^{n+1}-2^{n+1}=3\cdot3^n-2N.
\]

這還不含向量和、係數或索引 metadata。若按上述方式覆寫，所有階段可共用長度 **3ⁿ 的 work array**。4×4 的 block 層有效長度是 **9 → 6 → 4**，不需要同時保存三層的 19 個 byte。

這是此逐層排程的充分空間界，不是全域最小值證明。現有 C++ API 另配置 N bytes 輸出，因此每次額外 heap 是 3ⁿ+N；係數與 layout 另計。若改變全域計算順序，例如用顯式 stack 做深度優先排程，可以進一步探索 O(N) 的額外工作空間，但那已不是「先完成整層再到下一層」的排程，需另行實作與量測。

## 將論文 pseudo code 補成可執行規格

```text
layout = Compile(N, structure)
    建立唯一 V 的 ID 與兩個子和的相依關係
    建立唯一乘積及其係數 ID
    對 k=1..n 建立 block ID 與 F/G/H 子 ID
    決定拓樸順序、last-use 及可覆寫槽位

coefficients = PrepareMatrix(layout, A)
    計算共用係數子式與每個乘積需要的 X
    若啟用 unit-sum 特例，驗證條件

Evaluate(layout, coefficients, b)
    依順序計算 V；每個 ID 只算一次
    計算各唯一乘積，允許覆寫其 V
    對 k=1..n：
        每個唯一 block 只合併一次
        在最後一次讀取後回收／覆寫子值
    回傳 W[n,0,0]，長度 N
```

`PARAMS[-1].value` 因此可明確解讀為最後建立的 root block W[n,0,0] 的向量值。正式實作應保存 root ID，而非依賴容器剛好以這個順序排列。

## 已做的獨立驗證

執行 `python tools/derive_params.py`，以舊碼的 **lowbit** 生成順序重新建立公式，再核對本文件定義，涵蓋 N=1、2、4、8、16、32、64，三種矩陣情境，共 21 組。

N=1 是推導補上的合法 base case；舊 C++ 程式本身在 p=0 會越界，不能將此驗證解讀為舊程式已支援該邊界。

- 每個 leaf pos 的係數與向量索引集合一致。
- 同 pos 的對齊 block 具有完全相同的項排列；左上／左下共享成立。
- 父子的封閉索引式與三進位 ID 的 3t、3t+1、3t+2 關係成立。
- 每層唯一 block 數量、向量 XOR、輸出 XOR、一般與 unit 乘法數均符合上式。
- 將 root 展成 aᵢbⱼ 符號並消去重複項，恰好得到 brute-force 矩陣向量乘法。
- 獨立多項式乘法 oracle 的數值結果與現有 iterative 結果一致。

結果儲存在 [audit.json](../analysis/params-derivation/audit.json)。這是有限尺寸的符號與數值核對；一般性的理由是上述 pos 的高低位元分離、三子區塊恆等式與數學歸納。測試不是全域效能最優性的證明。

這次已將 `pos`、block 的 value、父子排程與回收規則連到現有迴圈。剩下屬於工程選擇的是係數子式去重的程度、metadata 的壓縮方式、是否採不同的全域排程；不用再把這些選擇誤當成公式尚未確定。
