#pragma once

#include <FL/Fl_Box.H>
#include <FL/Fl.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_ask.H>

#include <vector>

namespace LBFL
{
/**
 * A group that displays child controls with a margin around them.
 */
class MarginLayout : public Fl_Group
{
protected:
	/**
	 * Track the controls resident in the interior of this control.
	 * These are overlaid, so you can use transparency to build layered
	 * constructs.
	 */
	std::vector<Fl_Widget *> children;

	/** Left margin offset (pixels). */
	int leftMargin;

	/** Top margin offset (pixels). */
	int topMargin;

	/** Right margin offset (pixels). */
	int rightMargin;

	/** Bottom margin offset (pixels). */
	int bottomMargin;

public:
	/**
	 * Constructs a new BorderLayout, which emulates a Java Container managed
	 * by a java.awt.BorderLayout layout manager.
	 *
	 * <b>NB:</b> Unlike the FLTK groups, this container will not automatically
	 * snaffle newly constructed Widgets; the reason is that we cannot reliably
	 * <code>add()</code> widgets without specifying the location parameter.
	 *
	 * @param x The left edge of the control, nominally relative to the window.
	 * @param y The top edge of the control, nominally relative to the window.
	 * @param w The nominal width of the control.
	 * @param h The nominal height of the control.
	 * @param leftMargin The offset of the left edge of nested controls within
	 * 	this container.
	 * @param topMargin The offset of the top edge of nested controls within
	 * 	this container.
	 * @param rightMargin The offset of the right edge of nested controls within
	 * 	this container.
	 * @param bottomMargin The offset of the bottom edge of nested controls
	 * within this container.
	 */
	MarginLayout(int x, int y, int w, int h, int leftMargin = 5,
				 int topMargin = 5, int rightMargin = 5, int bottomMargin = 5)
		: Fl_Group(x, y, w, h),
		  leftMargin(leftMargin),
		  topMargin(topMargin),
		  rightMargin(rightMargin),
		  bottomMargin(bottomMargin)
	{
		end();
	}

	/**
	 *	<summary>Adds a child control at the centre of the container.
	 *	<list>
	 *	<item>This hides the <tt>Fl_Group::add</tt> method.</item>
	 *	</list>
	 *	</summary>
	 *	<param name="control">Reference to the new child control to add.</param>
	 */
	void add(Fl_Widget &control) /*override*/ { add(&control); }

	/**
	 *	<summary>Adds a child control at the centre of the container.
	 *	<list>
	 *	<item>This hides the <tt>Fl_Group::add</tt> method.</item>
	 *	</list>
	 *	</summary>
	 *	<param name="control">Address of the new child control to add.</param>
	 */

	void add(Fl_Widget *control) /*override*/
	{
		Fl_Group::add(control);
		children.insert(children.begin(), control);
		ResizeImpl();
		redraw();
	}

	/**
	 * <summary>
	 * 	Defines new location and dimensions for the container, and lays out the
	 * child controls using nominal dimensions obtained by querying their
	 * <c>w()</c> and <c>h()</c> properties. <list> <item>Overrides
	 * fl_Group::resize.</item> <item> Uses <c>ResizeImpl</c> to reposition and
	 * size each child based on its location and innate dimensions.
	 * 		</item>
	 * 	</list>
	 * </summary>
	 * <param name="xNew">The new left edge of the container.</param>
	 * <param name="yNew">The new top edge of the container.</param>
	 * <param name="wNew">The new width of the container.</param>
	 * <param name="hNew">The new height of the container.</param>
	 */

	void resize(int xNew, int yNew, int wNew, int hNew) override
	{
		// make new (x,y,w,h) values visible for children
		Fl_Widget::resize(xNew, yNew, wNew, hNew);
		ResizeImpl();
	}

	/**
	 *	<summary>
	 *	Lays out the child controls using nominal dimensions
	 *	obtained by querying their <c>w()</c> and <c>h()</c> properties.
	 *	</summary>
	 */
	virtual void ResizeImpl()
	{
		for (auto child : children)
		{
			child->resize(x() + leftMargin, y() + topMargin, w() - leftMargin - rightMargin, h() - topMargin - bottomMargin);
		}
	}

	void remove(Fl_Widget &child) { remove(&child); }

	void remove(Fl_Widget *child)
	{
		Fl_Group::remove(*child);

		auto iter = children.begin();

		while (iter != children.end())
		{
			if (*iter == child)
			{
				children.erase(iter);
				break;
			}
		}
	}
};
}  // namespace HBFL
