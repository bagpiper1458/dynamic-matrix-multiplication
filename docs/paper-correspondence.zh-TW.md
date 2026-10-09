# Brute force、Recursive、Iterative：演算法與論文對照

這份核對使用作者提供的 LaTeX 原稿，重點是「三子問題分解」、一般解的 `SM_A/SM_B`，以及 Algorithms 1–4。

**核對結論：recursive 的數學分解直接對應論文；目前 iterative 是相同一般解與共用計算概念的等價排程實作，但並非 Algorithms 1–4 的逐行翻譯。** 係數 DAG、三進位索引、連續 byte array、in-place 覆寫和跨矩陣共用 layout，都是實作時明確採用的工程選擇。作者後來確認舊碼僅作公式生成／驗證，PARAMS 共用運算當時只有 pseudo code。

後續的 [PARAMS 推導](params-derivation.zh-TW.md) 已將 pos、對齊 blocks 與三子合併的精確關係展開，提供現有迴圈排程與論文 block 概念之間的具體對照。

## 共同問題與記號

N=2ⁿ，計算 **N×N 矩陣乘 N×1 向量**。所有元素在 GF(2⁸)，加法是 XOR（⊕），乘法是多項式乘法模 0x11B。

Toeplitz 使用 A[i,j]=a[i−j]；circulant 再加上 a[k]=column[k mod N]。負對角索引對應 first row，非負索引對應 first column。輸入沒有展開成 N×N dense matrix，以免把配置 dense matrix 的額外成本算給 brute force。

## 1. Brute force：逐列逐項

對每個輸出直接計算：

\[
y_i=\bigoplus_{j=0}^{N-1}a_{i-j}\,b_j.
\]

```text
for i = 0 .. N-1:
    y[i] = mul(a[i], b[0])
    for j = 1 .. N-1:
        y[i] ^= mul(a[i-j], b[j])  # circulant 的係數索引取 mod N
```

共 N² 次 GF 乘法、N(N−1) 次 XOR；O(N²) 時間。C++ 輸出需要 N bytes，其他運算暫存是常數空間。它沒有預先建立的 layout 或矩陣 plan。

程式：`cpp/dynamic_matrix.hpp` 的 `toeplitz_matvec` / `circulant_matvec`。benchmark 中的 `baseline` 就是這個版本，圖表已統一標成 **Brute force**。

## 2. Recursive：論文的 F/G/H 三子問題

將 Toeplitz 矩陣與向量分成兩半：

\[
A=\begin{bmatrix}A_0&A_1\\A_2&A_0\end{bmatrix},\qquad
B=\begin{bmatrix}B_0\\B_1\end{bmatrix}.
\]

用 GF(2⁸) 中 X⊕X=0 的性質改寫：

\[
F=A_0(B_0\oplus B_1),\quad
G=(A_0\oplus A_1)B_1,\quad
H=(A_0\oplus A_2)B_0,
\]

\[
AB=\begin{bmatrix}F\oplus G\\F\oplus H\end{bmatrix}.
\]

例如上半部展開為 A₀B₀⊕A₀B₁⊕A₀B₁⊕A₁B₁，兩個 A₀B₁ 抵消，剩下原本的 A₀B₀⊕A₁B₁。下半部同理。A₀、A₀⊕A₁、A₀⊕A₂ 仍是 Toeplitz，因此同樣分解可持續套用，直到 1×1。

```text
evaluate(node, b):
    if size == 1:
        return mul(node.coefficient, b[0])
    split b into b0, b1
    f = evaluate(node.F, b0 XOR b1)
    g = evaluate(node.G, b1)
    h = evaluate(node.H, b0)
    return concat(f XOR g, f XOR h)
```

程式將「係數」與「向量運算」分開：`detail::compile` 先建立三子樹、保存葉係數；`detail::evaluate` 才對每個向量遞迴。保存整棵樹、預先算係數及共用 scratch 是實作選擇；F/G/H 的代數分解則直接來自論文。

設 M(n)、S(n) 是每次向量計算的 GF 乘法與 XOR 次數：

\[
M(n)=3M(n-1),\quad M(0)=1\quad\Rightarrow\quad M(n)=3^n,
\]

\[
S(n)=3S(n-1)+3\cdot2^{n-1},\quad S(0)=0
\quad\Rightarrow\quad S(n)=3(3^n-2^n).
\]

這些數字**不含係數預計算**。目前 recursive compile 另外執行 3·3ⁿ−4N+1 次係數 XOR。C++ 保存 (3·3ⁿ−1)/2 個節點；每次 apply 額外配置 N bytes 輸出及 2N−2 bytes scratch，還有 O(log N) 的呼叫堆疊。

對 circulant，最上層 A₁=A₂；遞迴過程仍用一般 Toeplitz 分解。如果列 XOR=1，沿每層只選 G/H 的 N 個葉係數為 1，可省下 N 次 GF 乘法，得到 3ⁿ−N。程式只根據這個結構條件省略乘法，不會把資料偶然等於 0 或 1 也當成一般性的省略。

## 3. Iterative：先共用一般解的項，再以迴圈計算

論文 Algorithm 1/2 先產生 `SM_A`、`SM_B`，將集合換成係數／向量元素的和，反轉 `SM_B` 的欄配對後得到輸出。Algorithm 3/4 再透過 `PARAMS(pos,value)` 共用已計算的 y 項與 blocks。

若真的把每列 N 個項全部重新計算，一共有 N² 個乘積出現位置。但其中有重複項。一般情況實際需要的不同乘積只有 **3ⁿ 個**；目前程式直接儲存這些不同乘積的排程。

### 一般解的精確表示

以下使用 0-based 輸出 i、欄遮罩 j。s⊆j 表示 s 的所有 1-bit 都是 j 的 1-bit；補集只取 n 個 bits：j̄=(N−1) XOR j。

\[
C(i,j)=\bigoplus_{s\subseteq j}a_{(i\mathbin{\&}j)-s},
\qquad
V(i,j)=\bigoplus_{t\subseteq\bar j}b_{((N-1-i)\mathbin{\&}j)\,|\,t}.
\]

\[
y_i=\bigoplus_{j=0}^{N-1} C(i,j)V(i,j).
\]

C 對應 `SM_A` 的集合和；V 對應**已反轉欄序**的 `SM_B` 集合和。circulant 將 a 的索引取 mod N。

一個 bit 只需記錄三種狀態：

| 狀態 | 此 bit 的係數選擇 | 向量索引 |
|---|---|---|
| j bit=0 | 不加入此位移 | 此 bit 自由取 0、1，兩者相加 |
| j bit=1、i bit=0 | 加入負位移 | 固定為 1 |
| j bit=1、i bit=1 | 加入正位移 | 固定為 0 |

程式把它們編為三進位數字 0、1、2；n 個 bits 共 3ⁿ 種狀態。這是壓縮一般解重複項的索引方法，論文沒有規定必須用這個資料結構。

### 程式的三階段

1. `IterativeLayout`：只由 N 與 Toeplitz/circulant 結構建立係數 XOR 關係、3ⁿ 個葉索引及 unit 葉位置。相同 operand-ID pair 只保留一次。這個 hash 去重只共用完全相同的二元子式，不保證找出所有代數等價式。
2. `layout.prepare`：將矩陣元素代入係數表，預計算 C 值並收集成長度 3ⁿ 的 byte array。這對應把矩陣相關的項預先算好。plan 之後不依賴 layout 的生命週期。
3. `plan.matvec`：逐層執行 `[u,v] → [u⊕v,v,u]`，建立不同 V 值；逐項乘上預計算係數；最後逐層執行 `[p₀,p₁,p₂] → [p₀⊕p₁,p₀⊕p₂]`，組合出 y。沒有 recursive calls。

向量展開用 3ⁿ−N 次 XOR，輸出合併用 2(3ⁿ−N) 次 XOR，合計仍為 3(3ⁿ−N)。同一個乘積或合併結果直接供多個輸出共用，不逐列重算。這保留了 PARAMS 避免重算的概念；固定位置的陣列與迴圈取代了執行時搜尋 `PARAMS(pos,value)`。

### 4×4 完整乘積對照

令 v01=b₀⊕b₁、v23=b₂⊕b₃、vAll=v01⊕v23、v13=b₁⊕b₃、v02=b₀⊕b₂。這五個 XOR 只算一次。

| 乘積 | 預計算的係數 | 每次的向量項 |
|---|---|---|
| p₀ | a₀ | vAll |
| p₁ | a₀⊕a₋₁ | v13 |
| p₂ | a₀⊕a₁ | v02 |
| p₃ | a₀⊕a₋₂ | v23 |
| p₄ | a₀⊕a₋₁⊕a₋₂⊕a₋₃ | b₃ |
| p₅ | a₀⊕a₁⊕a₋₂⊕a₋₁ | b₂ |
| p₆ | a₀⊕a₂ | v01 |
| p₇ | a₀⊕a₋₁⊕a₂⊕a₁ | b₁ |
| p₈ | a₀⊕a₁⊕a₂⊕a₃ | b₀ |

共用合併項：

```text
f0 = p0 XOR p1     f1 = p0 XOR p2
g0 = p3 XOR p4     g1 = p3 XOR p5
h0 = p6 XOR p7     h1 = p6 XOR p8
y0 = f0 XOR g0    y1 = f1 XOR g1
y2 = f0 XOR h0    y3 = f1 XOR h1
```

所以是 9 次乘法、5+10=15 次向量端 XOR。unit-sum circulant 的 p₄、p₅、p₇、p₈ 係數為 1，剩 5 次乘法。4×4 係數預計算則為 Toeplitz 11 XOR、circulant 5 XOR；這些是目前去重排程的數字，不能混入論文每個向量的 15 XOR。

## 哪些相同，哪些不是逐行重現

| 論文內容 | 目前程式對應 | 核對狀態 |
|---|---|---|
| F/G/H 三子問題 | `detail::compile/evaluate` | 直接使用相同分解 |
| Algorithm 1/2 的 SM 集合與欄反轉 | 上述 C/V 一般式、三進位 product IDs | 公式等價；主實作不建立 N×N SM 表 |
| PARAMS 共用 y 項與 block 結果 | 向量展開、唯一乘積、共用合併結果 | 相同共用目的；陣列位置與遍歷順序由重建指定 |
| 矩陣固定時可預計算 | `layout.prepare` | 明確分離；額外量測準備成本 |
| 係數子式去重 | canonical operand pair 的 hash 表 | 工程選擇；原稿沒有給出此 hash 演算法 |
| 暫存空間生命週期 | 長度 3ⁿ 的原地 work array | 工程選擇，未證明最小 |
| GF(2⁸) scalar multiplication | `gf_mul`，模 0x11B | 對應原稿提供的 scalar code |

原稿對 `SM_B` 初值在文字敘述中出現 bᵢ 與 b[N−1−i] 兩種寫法；此處採用 Algorithms 1/2 所列的 b[N−1−i]，並配合欄反轉，與 4×4 展開一致。Algorithm 4 的產生器呼叫寫為 Algorithm1，而前文指定 circulant 的 Algorithm2；目前按前文加上負索引 mod N 處理。這些採用方式公開記錄，不能當成原稿每一行都無歧義。

## 如何驗證，而不是只比較最終答案

- `tests/test_paper_correspondence.py` 獨立依 Algorithm 1/2 的「前集合＋位移集合」生成 SM，使用集合運算，不沿用主程式的三進位生成器。
- N=1、2、4、8、16、32：將主程式每個係數 XOR 展成符號，逐項對上 SM；再將每列的乘積展開成 a[k]·b[j]，以 GF(2) 消去重複項，結果必須恰好是 brute force 的那一列。這是這些尺寸的符號核對，涵蓋任意 byte 輸入；不是單一隨機值比對，也不是對所有 N 的機器證明。
- Python 與 C++ 都另以 Algorithm 1/2 的獨立數值計算，核對 N≤16 的 recursive 與 iterative 路徑。C++ 額外核對已確認的 4×4 九個係數順序。
- 既有 C++ 三版本測試覆蓋 N=1 到 1024 的二次方尺寸，另有 72 組共用測資、65,536 個有限體乘積的獨立 polynomial oracle。

代數上的一般性依據仍是 F/G/H 恆等式及每一位的三種狀態分解。測試提供實作證據，不用測試範圍取代數學證明，也不把輸出一致誇大為 pseudo code 已唯一指定儲存排程。

## 三版本比較的正確口徑

| 每個向量，N=2ⁿ | Brute force | Recursive | Iterative |
|---|---|---|---|
| 一般乘法次數 | N² | 3ⁿ | 3ⁿ |
| unit-sum circulant 乘法次數 | N² | 3ⁿ−N | 3ⁿ−N |
| 向量端 XOR | N(N−1) | 3(3ⁿ−N) | 3(3ⁿ−N) |
| 已準備後的時間階 | O(N²) | O(N^log₂3) | O(N^log₂3) |
| 每次額外 heap，含輸出 | N bytes | 3N−2 bytes | 3ⁿ+N bytes |

最後一列不含保存的係數、layout、輸入、stack 或 allocator overhead；它不能取代總記憶體比較。首次建立、只共用尺寸排程、共用矩陣三種情境的完整時間與記憶體結果，見 [三版本報告](../analysis/iterative/report.html)。
