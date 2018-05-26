# dom-sagecell

## Purpose

This is a repository to store python code for use in the sagecell server.

My python code is about euclidean geometry projects.

## Example of referencing python code in sagemath cell

Add a file *sample1.sage* in main directory of the github repository containing :

```python
"""
def f(x):
  return x**2
  a=3
"""
```

Click _Raw_ button, to display file content in new browser windows and copy URL.

Reference it in the cell of the sagecell server. Type in the cell :

```python
import urllib2
url = 'https://raw.githubusercontent.com/dominiquelaurain/dom-sagecell/master/sample1.sage'
exec(eval(urllib2.urlopen(url, None).read()))

print "a = ",a,"and f(a) = a^2 = ",f(a)
```


Then click _Evaluate_ button. You will get :

```text
a =  3 and f(a) = a^2 =  9
```

Global variable a and function f have been loaded into the cell after exec call.


## References

- Markdown - cheat sheet - https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet
- Sagecell server - sagemath cell - https://sagecell.sagemath.org/
- Google group sagecell - how to reference external code - https://groups.google.com/forum/#!topic/sage-cell/__tQ_iG9FZQ
