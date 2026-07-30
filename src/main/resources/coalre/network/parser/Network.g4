grammar Network;

network: node ';'? EOF;

node: ('(' node (',' node)* ')')? post ;

post: label? hybrid? meta? (':' length=number)? ;

label: number | string ;

hybrid: '#' (type=STRINGALPHA)? id=INT ;

meta: '[&' attrib (',' attrib)* ']' ;

attrib: attribKey=string '=' attribValue ;

attribValue: number | string | vector ;

number: INT | FLOAT ;

vector: '{' attribValue (',' attribValue)* '}' ;

string:  STRINGALPHA (STRINGNUMALPHA | INT)? | STRINGNUMALPHA;

// Lexer rules:

FLOAT : '-'? NNINT ('.' D*) ([eE] ('-'|'+')? D+)? ;
INT : '-'? NNINT;
fragment NNINT : '0' | NZD D* ;
fragment NZD : [1-9] ;
fragment D : [0-9] ;

STRINGALPHA: [a-zA-Z]+ ;

STRINGNUMALPHA :
    [0-9|*%/.\-+_&][a-zA-Z0-9|*%/.\-+_&]+  // these chars don't need quotes
    | '"' .*? '"'
    | '\'' .*? '\''
    ;

WHITESPACE : [ \t\r\n]+ -> skip ;