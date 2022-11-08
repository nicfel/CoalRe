// Generated from /Users/vaughant/code/beast_and_friends/CoalRe/src/coalre/network/parser/Network.g4 by ANTLR 4.10.1
package coalre.network.parser;
import org.antlr.v4.runtime.tree.ParseTreeListener;

/**
 * This interface defines a complete listener for a parse tree produced by
 * {@link NetworkParser}.
 */
public interface NetworkListener extends ParseTreeListener {
	/**
	 * Enter a parse tree produced by {@link NetworkParser#network}.
	 * @param ctx the parse tree
	 */
	void enterNetwork(NetworkParser.NetworkContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#network}.
	 * @param ctx the parse tree
	 */
	void exitNetwork(NetworkParser.NetworkContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#node}.
	 * @param ctx the parse tree
	 */
	void enterNode(NetworkParser.NodeContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#node}.
	 * @param ctx the parse tree
	 */
	void exitNode(NetworkParser.NodeContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#post}.
	 * @param ctx the parse tree
	 */
	void enterPost(NetworkParser.PostContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#post}.
	 * @param ctx the parse tree
	 */
	void exitPost(NetworkParser.PostContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#label}.
	 * @param ctx the parse tree
	 */
	void enterLabel(NetworkParser.LabelContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#label}.
	 * @param ctx the parse tree
	 */
	void exitLabel(NetworkParser.LabelContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#hybrid}.
	 * @param ctx the parse tree
	 */
	void enterHybrid(NetworkParser.HybridContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#hybrid}.
	 * @param ctx the parse tree
	 */
	void exitHybrid(NetworkParser.HybridContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#meta}.
	 * @param ctx the parse tree
	 */
	void enterMeta(NetworkParser.MetaContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#meta}.
	 * @param ctx the parse tree
	 */
	void exitMeta(NetworkParser.MetaContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#attrib}.
	 * @param ctx the parse tree
	 */
	void enterAttrib(NetworkParser.AttribContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#attrib}.
	 * @param ctx the parse tree
	 */
	void exitAttrib(NetworkParser.AttribContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#attribValue}.
	 * @param ctx the parse tree
	 */
	void enterAttribValue(NetworkParser.AttribValueContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#attribValue}.
	 * @param ctx the parse tree
	 */
	void exitAttribValue(NetworkParser.AttribValueContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#number}.
	 * @param ctx the parse tree
	 */
	void enterNumber(NetworkParser.NumberContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#number}.
	 * @param ctx the parse tree
	 */
	void exitNumber(NetworkParser.NumberContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#vector}.
	 * @param ctx the parse tree
	 */
	void enterVector(NetworkParser.VectorContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#vector}.
	 * @param ctx the parse tree
	 */
	void exitVector(NetworkParser.VectorContext ctx);
	/**
	 * Enter a parse tree produced by {@link NetworkParser#string}.
	 * @param ctx the parse tree
	 */
	void enterString(NetworkParser.StringContext ctx);
	/**
	 * Exit a parse tree produced by {@link NetworkParser#string}.
	 * @param ctx the parse tree
	 */
	void exitString(NetworkParser.StringContext ctx);
}