/*******************************************************************************
 * Copyright (c) 2016-2019 Embedded Systems and Applications Group
 * Department of Computer Science, Technische Universitaet Darmstadt,
 * Hochschulstr. 10, 64289 Darmstadt, Germany.
 *
 * All rights reserved.
 *
 * This software is provided free for educational use only.
 * It may not be used for commercial purposes without the
 * prior written permission of the authors.
 ******************************************************************************/
package mavlc.context_analysis;

import mavlc.errors.NonConstantExpressionError;
import mavlc.syntax.AstNode;
import mavlc.syntax.AstNodeBaseVisitor;
import mavlc.syntax.expression.*;

/* TODO enter group information
*
* EiCB group number: 20
* Names and matriculation numbers of all group members:
* Devran Akdag  2925536 
* Jehad Sefi    2430157
* Minh Huy Tran 2712710

*/

public class ConstantExpressionEvaluator extends AstNodeBaseVisitor<Integer, Void> {
	@Override
	protected Integer defaultOperation(AstNode node, Void obj) {
		if(node instanceof Expression) {
			throw new NonConstantExpressionError((Expression) node);
		} else {
			throw new RuntimeException("Internal compiler error: should not try to constant-evaluate non-expressions");
		}
	}
	
	@Override
	public Integer visitIntValue(IntValue intValue, Void __) {
		return intValue.value;
	}
	
	// TODO implement (exercise 2.3)
	@Override 
	public Integer visitAddition(Addition addition, Void __) { 
		Integer leftOperand = addition.leftOperand.accept(this);
		Integer rightOperand = addition.rightOperand.accept(this);
		return leftOperand + rightOperand;
	}
	
	@Override 
	public Integer visitSubtraction(Subtraction subtraction, Void __) { 
		Integer leftOperand = subtraction.leftOperand.accept(this);
		Integer rightOperand = subtraction.rightOperand.accept(this);
		return leftOperand - rightOperand;
	}
	
	@Override 
	public Integer visitMultiplication(Multiplication multiplication, Void __) { 
		Integer leftOperand = multiplication.leftOperand.accept(this);
		Integer rightOperand = multiplication.rightOperand.accept(this);
		return leftOperand * rightOperand;
	}
	
	@Override 
	public Integer visitDivision(Division division, Void __) { 
		Integer leftOperand = division.leftOperand.accept(this);
		Integer rightOperand = division.rightOperand.accept(this);
		return leftOperand / rightOperand;
	}
	
	@Override 
	public Integer visitExponentiation(Exponentiation exponentiation, Void __) { 
		Integer leftOperand = exponentiation.leftOperand.accept(this);
		Integer rightOperand = exponentiation.rightOperand.accept(this);
		return (int)Math.pow(leftOperand, rightOperand);
	}
	
	@Override 
	public Integer visitUnaryMinus(UnaryMinus unaryMinus, Void __) { 
		Integer operand =  unaryMinus.operand.accept(this);
		return -operand;
	}
}
