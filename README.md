***
## Gitの使い方
***

### GitHubで100 MB以上のファイルをpushしようとしてエラーになった後の解決方法
> 参考資料: <br>
> 　[aa](https://qiita.com/ffggss/items/b61e726bc8cbd3137956) <br>
> ログを見て、origin/master, origin/HEAD記載箇所のmd5値を指定して、コマンド実行 <br>
> ```bash
> $ git reset --soft
> ```

### Gitのレポジトリ管理が変更（2021/08/14）
> git push -f origin main でアップデートしようしたら、エラーメッセージが出現 <br>
> 解決策: <br>
> 個人アクセストークンを作成して、git push -f origin main実行時のパスワードに作成したトークンを入力した。 <br>
>
> 参考資料: <br>
> ・GitHubがパスワード認証を廃止するらしいので <br>
>　[https://qiita.com/shiro01/items/e886aa1e4beb404f9038]() <br>
> ・GitHubでhttpsのパスワード認証が廃止された。Please use a personal access token instead. <br>
>　[https://qiita.com/shunsa10/items/e43564cf48f84b95455b]() <br>
>
> Username: touchiyama <br>
> Personal access tokens(パスワード): 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 <br>
>
> Webで自分のアカウントに入る https://github.com/touchiyama　<br>
> レポジトリを作成する <br>

### 基本操作
> 以下、ターミナル上で操作 <br>
> ```bash
> $ cd /Volumes/AutoReg
> $ echo "# AutoReg" >> README.md
> $ git init
> $ git add README.md
> $ git commit -m "first commit" #必ずcommitさせる
> $ git branch -M main  #branchがmainになる
> $ git status #changed to commitならば、pushが可能な状態
> $ git remote add origin https://github.com/touchiyama/AutoReg.git
> $ git push -u origin main　#README.mdがAutoRegという名前のgitにpushされる
>
> $ git add .  # ’.’はカレントディレクトリにある全ての要素を示す
> $ git status　
> $ git commit -m "first commit"  #必ずcommitさせる　　　
> $ git push -f origin main  #最新バージョンにアップデートさせる
> ```
> 参考資料: <br>
>　[https://qiita.com/kazuki0714/items/ceda3a6721a9a99082de]() <br>