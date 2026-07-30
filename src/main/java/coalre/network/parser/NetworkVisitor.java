// Generated from java-escape by ANTLR 4.11.1
package coalre.network.parser;
import org.antlr.v4.runtime.tree.ParseTreeVisitor;

/**
 * This interface defines a complete generic visitor for a parse tree produced
 * by {@link NetworkParser}.
 *
 * @param <T> The return type of the visit operation. Use {@link Void} for
 * operations with no return type.
 */
public interface NetworkVisitor<T> extends ParseTreeVisitor<T> {
	/**
	 * Visit a parse tree produced by {@link NetworkParser#network}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitNetwork(NetworkParser.NetworkContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#node}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitNode(NetworkParser.NodeContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#post}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitPost(NetworkParser.PostContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#label}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitLabel(NetworkParser.LabelContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#hybrid}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitHybrid(NetworkParser.HybridContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#meta}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitMeta(NetworkParser.MetaContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#attrib}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitAttrib(NetworkParser.AttribContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#attribValue}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitAttribValue(NetworkParser.AttribValueContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#number}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitNumber(NetworkParser.NumberContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#vector}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitVector(NetworkParser.VectorContext ctx);
	/**
	 * Visit a parse tree produced by {@link NetworkParser#string}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitString(NetworkParser.StringContext ctx);
}