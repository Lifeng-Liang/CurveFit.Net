
CurveFit.Net
==========

This is a pure C# curve fit library. I translate it from java version [CurveFitter](https://github.com/BinRoot/CurveFitter).

这是一个纯 .Net 版的曲线拟合库，我从java版的 [CurveFitter](https://github.com/BinRoot/CurveFitter) 转换而来。

转换中，几乎没有修改任何东西，只是调整一些必要的地方，比如系统函数的大小写，关键字冲突，java特有函数之类的。

简单看了一下代码结构，其实不算好，不过测试了一下，用它计算出的幂函数拟合结果，和使用 python scipy 计算出的结果，非常接近，可以认为可用性还是非常高的了吧。

因为没有实际修改任何地方，所以也没有进行详细测试，基本可以认为，java版能做的，.net 版的表现都是一致的。
